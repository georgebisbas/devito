from __future__ import absolute_import

from sympy import Indexed

from devito.dimension import LoweredDimension
from devito.dle import retrieve_iteration_tree
from devito.logger import yask as log, yask_warning as warning
from devito.operator import OperatorRunnable
from devito.visitors import FindSymbols

from devito.yask.kernel import (cfac, nfac, ofac, namespace, yask_context,
                                yask_jit, _force_exit)

__all__ = ['Operator']


class Operator(OperatorRunnable):

    """
    A special :class:`OperatorCore` to JIT-compile and run operators through YASK.
    """

    def __init__(self, expressions, **kwargs):
        kwargs['dle'] = 'noop'
        super(Operator, self).__init__(expressions, **kwargs)

    def _specialize(self, nodes, elemental_functions):
        """
        Create a YASK representation of this Iteration/Expression tree.
        """

        log("Specializing a Devito Operator for YASK...")

        # Set up the YASK solution
        self.soln = cfac.new_solution(namespace['kernel-real'])

        # Silence YASK
        self.yask_compiler_output = ofac.new_string_output()
        self.soln.set_debug_output(self.yask_compiler_output)

        trees = retrieve_iteration_tree(nodes)
        if len(trees) > 1:
            _force_exit("Currently unable to handle Operators w/ more than 1 loop nests")

        tree = trees[0]
        candidate = tree[-1]

        expressions = [e for e in candidate.nodes if e.is_Expression]
        functions = [i for i in FindSymbols().visit(candidate) if i.is_TimeData]
        keys = set([(i.indices, i.shape, i.dtype, i.space_order) for i in functions])
        if len(keys) > 1:
            _force_exit("Currently unable to handle Operators w/ heterogeneous grids")
        self.context = yask_context(*keys.pop())

        # Perform the translation on an expression basis
        transform = sympy2yask(self.context, self.soln)
        try:
            for i in expressions:
                ast = transform(i.expr)
                log("Converted %s into YASK AST [%s]", str(i.expr), ast.format_simple())
        except:
            _force_exit("Couldn't convert %s into YASK format" % str(i.expr))

        # Print some useful information about the newly constructed Yask solution
        log("Solution '" + self.soln.get_name() + "' contains " +
            str(self.soln.get_num_grids()) + " grid(s), and " +
            str(self.soln.get_num_equations()) + " equation(s).")

        # Set necessary run-time parameters
        self.soln.set_step_dim_name(self.context.time_dimension)
        self.soln.set_domain_dim_names(self.context.space_dimensions)
        self.soln.set_element_bytes(4)

        # JIT YASK kernel
        name = 'kernel%d' % self.context.nkernels
        yk = yask_jit(self.soln, name)

        # Build an empty kernel solution
        kfac = yk.yk_factory()
        env = kfac.new_env()
        self.ksoln = kfac.new_solution(env)

        # Silence YASK
        self.yask_kernel_output = yk.yask_output_factory().new_string_output()
        self.ksoln.set_debug_output(self.yask_kernel_output)

        self.ksoln.set_num_ranks(self.ksoln.get_domain_dim_name(0), env.get_num_ranks())

        # Track this Operator
        self.context.add_kernel(name, self)

        # TODO: need to update nodes and elemental_functions

        log("Specialization successfully performed!")

        return nodes, elemental_functions

    def apply(self, *args, **kwargs):
        # Build the arguments list to invoke the kernel function
        arguments, dim_sizes = self.arguments(*args, **kwargs)

        # Set the domain sizes
        for dm in self.ksoln.get_domain_dim_names():
            self.ksoln.set_rank_domain_size(dm, dim_sizes[dm])

        # Share the grids from the hook solution
        for kgrid in self.ksoln.get_grids():
            name = kgrid.get_name()
            try:
                hook_grid = self.context.grids[name]
            except KeyError:
                _force_exit("Unknown grid %s" % name)
            kgrid.share_storage(hook_grid)
            log("Shared storage from grid <%s>" % hook_grid.get_name())

        # Print some info about the solution.
        log("YASK Stencil-solution '" + self.ksoln.get_name() + "':")
        log("  Step dimension: " + repr(self.ksoln.get_step_dim_name()))
        log("  Domain dimensions: " + repr(self.ksoln.get_domain_dim_names()))
        log("  Grids:")
        for grid in self.ksoln.get_grids():
            log("    " + grid.get_name() + repr(grid.get_dim_names()))

        log("Running Operator through YASK...")
        self.ksoln.prepare_solution()
        # TODO: getting number of timesteps in a hacky way
        self.ksoln.run_solution(arguments["%s_size" % self.context.time_dimension])


class sympy2yask(object):
    """
    Convert a SymPy expression into a YASK abstract syntax tree and create any
    necessay YASK grids.
    """

    def __init__(self, context, soln):
        self.context = context
        self.soln = soln
        self.mapper = {}

    def __call__(self, expr):

        def nary2binary(args, op):
            r = run(args[0])
            return r if len(args) == 1 else op(r, nary2binary(args[1:], op))

        def run(expr):
            if expr.is_Integer:
                return nfac.new_const_number_node(int(expr))
            elif expr.is_Float:
                return nfac.new_const_number_node(float(expr))
            elif expr.is_Symbol:
                assert expr in self.mapper
                return self.mapper[expr]
            elif isinstance(expr, Indexed):
                function = expr.base.function
                name = function.name
                if name not in self.context.grids:
                    function.data  # Create uninitialized grid (i.e., 0.0 everywhere)
                if name not in self.mapper:
                    dimensions = self.context.grids[name].get_dim_names()
                    self.mapper[name] = self.soln.new_grid(name, dimensions)
                indices = [int((i.origin if isinstance(i, LoweredDimension) else i) - j)
                           for i, j in zip(expr.indices, function.indices)]
                return self.mapper[name].new_relative_grid_point(*indices)
            elif expr.is_Add:
                return nary2binary(expr.args, nfac.new_add_node)
            elif expr.is_Mul:
                return nary2binary(expr.args, nfac.new_multiply_node)
            elif expr.is_Pow:
                num, den = expr.as_numer_denom()
                if num == 1:
                    return nfac.new_divide_node(run(num), run(den))
            elif expr.is_Equality:
                if expr.lhs.is_Symbol:
                    assert expr.lhs not in self.mapper
                    self.mapper[expr.lhs] = run(expr.rhs)
                else:
                    return nfac.new_equation_node(*[run(i) for i in expr.args])
            else:
                warning("Missing handler in Devito-YASK translation")
                raise NotImplementedError

        return run(expr)
