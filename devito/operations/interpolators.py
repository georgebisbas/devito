from abc import ABC, abstractmethod

import sympy
import numpy as np
from cached_property import cached_property

from devito.symbolics import retrieve_function_carriers
from devito.tools import as_tuple, powerset, flatten, prod
from devito.types import (ConditionalDimension, Eq, Inc, Evaluable, Symbol)

__all__ = ['LinearInterpolator', 'PrecomputedInterpolator']


class UnevaluatedSparseOperation(sympy.Expr, Evaluable):

    """
    Represents an Injection or an Interpolation operation performed on a
    SparseFunction. Evaluates to a list of Eq objects.

    Parameters
    ----------
    interpolator : Interpolator
        Interpolator object that will be used to evaluate the operation.
    callback : callable
        A routine generating the symbolic expressions for the operation.
    """

    subdomain = None

    def __new__(cls, interpolator, callback):
        obj = super().__new__(cls)

        obj.interpolator = interpolator
        obj.callback = callback

        return obj

    def _evaluate(self, **kwargs):
        return_value = self.callback()
        assert(all(isinstance(i, Eq) for i in return_value))
        return return_value

    def __add__(self, other):
        return flatten([self, other])

    def __radd__(self, other):
        return flatten([other, self])


class Interpolation(UnevaluatedSparseOperation):

    """
    Represents an Interpolation operation performed on a SparseFunction.
    Evaluates to a list of Eq objects.
    """

    def __new__(cls, expr, offset, increment, self_subs, interpolator, callback):
        obj = super().__new__(cls, interpolator, callback)

        # TODO: unused now, but will be necessary to compute the adjoint
        obj.expr = expr
        obj.offset = offset
        obj.increment = increment
        obj.self_subs = self_subs

        return obj

    def __repr__(self):
        return "Interpolation(%s into %s)" % (repr(self.expr),
                                              repr(self.interpolator.sfunction))


class Injection(UnevaluatedSparseOperation):

    """
    Represents an Injection operation performed on a SparseFunction.
    Evaluates to a list of Eq objects.
    """

    def __new__(cls, field, expr, offset, interpolator, callback):
        obj = super().__new__(cls, interpolator, callback)

        # TODO: unused now, but will be necessary to compute the adjoint
        obj.field = field
        obj.expr = expr
        obj.offset = offset

        return obj

    def __repr__(self):
        return "Injection(%s into %s)" % (repr(self.expr), repr(self.field))


class GenericInterpolator(ABC):

    """
    Abstract base class defining the interface for an interpolator.
    """

    @abstractmethod
    def inject(self, *args, **kwargs):
        pass

    @abstractmethod
    def interpolate(self, *args, **kwargs):
        pass


class LinearInterpolator(GenericInterpolator):

    """
    Concrete implementation of GenericInterpolator implementing a Linear interpolation
    scheme, i.e. Bilinear for 2D and Trilinear for 3D problems.

    Parameters
    ----------
    sfunction: The SparseFunction that this Interpolator operates on.
    """

    def __init__(self, sfunction):
        self.sfunction = sfunction

    @property
    def grid(self):
        return self.sfunction.grid

    @cached_property
    def _interpolation_coeffs(self):
        """
        Symbolic expression for the coefficients for sparse point interpolation
        according to:

            https://en.wikipedia.org/wiki/Bilinear_interpolation.

        Returns
        -------
        Matrix of coefficient expressions.
        """
        # Grid indices corresponding to the corners of the cell ie x1, y1, z1
        indices1 = tuple(sympy.symbols('%s1' % d) for d in self.grid.dimensions)
        indices2 = tuple(sympy.symbols('%s2' % d) for d in self.grid.dimensions)
        # 1, x1, y1, z1, x1*y1, ...
        indices = list(powerset(indices1))
        indices[0] = (1,)
        point_sym = list(powerset(self.sfunction._point_symbols))
        point_sym[0] = (1,)
        # 1, px. py, pz, px*py, ...
        A = []
        ref_A = [np.prod(ind) for ind in indices]
        # Create the matrix with the same increment order as the point increment
        for i in self.sfunction._point_increments:
            # substitute x1 by x2 if increment in that dimension
            subs = dict((indices1[d], indices2[d] if i[d] == 1 else indices1[d])
                        for d in range(len(i)))
            A += [[1] + [a.subs(subs) for a in ref_A[1:]]]

        A = sympy.Matrix(A)
        # Coordinate values of the sparse point
        p = sympy.Matrix([[np.prod(ind)] for ind in point_sym])

        # reference cell x1:0, x2:h_x
        left = dict((a, 0) for a in indices1)
        right = dict((b, dim.spacing) for b, dim in zip(indices2, self.grid.dimensions))
        reference_cell = {**left, **right}
        # Substitute in interpolation matrix
        A = A.subs(reference_cell)
        return A.inv().T * p

    def _interpolation_indices(self, variables, offset=0, field_offset=0,
                               implicit_dims=None):
        """
        Generate interpolation indices for the DiscreteFunctions in ``variables``.
        """
        index_matrix, points = self.sfunction._index_matrix(offset)

        idx_subs = []
        for i, idx in enumerate(index_matrix):
            # Introduce ConditionalDimension so that we don't go OOB
            mapper = {}
            for j, d in zip(idx, self.grid.dimensions):
                p = points[j]
                lb = sympy.And(p >= d.symbolic_min - self.sfunction._radius,
                               evaluate=False)
                ub = sympy.And(p <= d.symbolic_max + self.sfunction._radius,
                               evaluate=False)
                condition = sympy.And(lb, ub, evaluate=False)
                mapper[d] = ConditionalDimension(p.name, self.sfunction._sparse_dim,
                                                 condition=condition, indirect=True)

            # Apply mapper to each variable with origin correction before the
            # Dimensions get replaced
            subs = {v: v.subs({k: c - v.origin.get(k, 0) for k, c in mapper.items()})
                    for v in variables}

            # Track Indexed substitutions
            idx_subs.append(subs)

        # Temporaries for the position
        temps = [Eq(v, k, implicit_dims=implicit_dims)
                 for k, v in self.sfunction._position_map.items()]
        # Temporaries for the indirection dimensions
        temps.extend([Eq(v, k.subs(self.sfunction._position_map),
                         implicit_dims=implicit_dims)
                      for k, v in points.items()])
        # Temporaries for the coefficients
        temps.extend([Eq(p, c.subs(self.sfunction._position_map),
                         implicit_dims=implicit_dims)
                      for p, c in zip(self.sfunction._point_symbols,
                                      self.sfunction._coordinate_bases(field_offset))])

        return idx_subs, temps

    def subs_coords(self, _expr, *idx_subs):
        return [_expr.xreplace(v_sub) * b.xreplace(v_sub)
                for b, v_sub in zip(self._interpolation_coeffs, idx_subs)]

    def subs_coords_eq(self, field, _expr, *idx_subs, implicit_dims=None):
        return [Inc(field.xreplace(vsub), _expr.xreplace(vsub) * b,
                    implicit_dims=implicit_dims)
                for b, vsub in zip(self._interpolation_coeffs, idx_subs)]

    def implicit_dims(self, implicit_dims):
        return as_tuple(implicit_dims) + self.sfunction.dimensions

    def interpolate(self, expr, offset=0, increment=False, self_subs={},
                    implicit_dims=None):
        """
        Generate equations interpolating an arbitrary expression into ``self``.

        Parameters
        ----------
        expr : expr-like
            Input expression to interpolate.
        offset : int, optional
            Additional offset from the boundary.
        increment: bool, optional
            If True, generate increments (Inc) rather than assignments (Eq).
        implicit_dims : Dimension or list of Dimension, optional
            An ordered list of Dimensions that do not explicitly appear in the
            interpolation expression, but that should be honored when constructing
            the operator.
        """
        implicit_dims = self.implicit_dims(implicit_dims)

        def callback():
            # Derivatives must be evaluated before the introduction of indirect accesses
            try:
                _expr = expr.evaluate
            except AttributeError:
                # E.g., a generic SymPy expression or a number
                _expr = expr

            variables = list(retrieve_function_carriers(_expr))

            # Need to get origin of the field in case it is staggered
            # TODO: handle each variable staggereing spearately
            field_offset = variables[0].origin
            # List of indirection indices for all adjacent grid points
            idx_subs, temps = self._interpolation_indices(
                variables, offset, field_offset=field_offset, implicit_dims=implicit_dims
            )

            # Substitute coordinate base symbols into the interpolation coefficients
            args = self.subs_coords(_expr, *idx_subs)

            # Accumulate point-wise contributions into a temporary
            rhs = Symbol(name='sum', dtype=self.sfunction.dtype)
            summands = [Eq(rhs, 0., implicit_dims=implicit_dims)]
            summands.extend([Inc(rhs, i, implicit_dims=implicit_dims) for i in args])

            # Write/Incr `self`
            lhs = self.sfunction.subs(self_subs)
            ecls = Inc if increment else Eq
            last = [ecls(lhs, rhs, implicit_dims=implicit_dims)]

            return temps + summands + last

        return Interpolation(expr, offset, increment, self_subs, self, callback)

    def inject(self, field, expr, offset=0, implicit_dims=None):
        """
        Generate equations injecting an arbitrary expression into a field.

        Parameters
        ----------
        field : Function
            Input field into which the injection is performed.
        expr : expr-like
            Injected expression.
        offset : int, optional
            Additional offset from the boundary.
        implicit_dims : Dimension or list of Dimension, optional
            An ordered list of Dimensions that do not explicitly appear in the
            injection expression, but that should be honored when constructing
            the operator.
        """
        implicit_dims = self.implicit_dims(implicit_dims)

        def callback():
            # Derivatives must be evaluated before the introduction of indirect accesses
            try:
                _expr = expr.evaluate
            except AttributeError:
                # E.g., a generic SymPy expression or a number
                _expr = expr

            variables = list(retrieve_function_carriers(_expr)) + [field]

            # Need to get origin of the field in case it is staggered
            field_offset = field.origin
            # List of indirection indices for all adjacent grid points
            idx_subs, temps = self._interpolation_indices(
                variables, offset, field_offset=field_offset, implicit_dims=implicit_dims
            )

            # Substitute coordinate base symbols into the interpolation coefficients
            eqns = self.subs_coords_eq(field, _expr, *idx_subs,
                                       implicit_dims=implicit_dims)

            return temps + eqns

        return Injection(field, expr, offset, self, callback)


class PrecomputedInterpolator(LinearInterpolator):

    def __init__(self, obj):
        self.sfunction = obj

    @property
    def r(self):
        return self.obj.r

    def _interpolation_indices(self, variables, offset=0, field_offset=0,
                               implicit_dims=None):
        """
        Generate interpolation indices for the DiscreteFunctions in ``variables``.
        """
        if self.sfunction.gridpoints is None:
            return super()._interpolation_indices(variables, offset=offset,
                                                  field_offset=field_offset,
                                                  implicit_dims=implicit_dims)

        index_matrix, points, shifts = self.sfunction._index_matrix(offset)

        idx_subs = []
        coeffs = self._interpolation_coeffs
        dt, it = coeffs.dimensions[1:]
        for i, idx in enumerate(index_matrix):
            # Introduce ConditionalDimension so that we don't go OOB
            mapper = {}
            for j, (di, d) in zip(idx, enumerate(self.grid.dimensions)):
                p = points[j]
                lb = sympy.And(p >= d.symbolic_min - self.sfunction.r // 2,
                               evaluate=False)
                ub = sympy.And(p <= d.symbolic_max + self.sfunction.r // 2,
                               evaluate=False)
                condition = sympy.And(lb, ub, evaluate=False)
                mapper[d] = ConditionalDimension(p.name, self.sfunction._sparse_dim,
                                                 condition=condition, indirect=True)
                mapper[coeffs._subs(dt, di)] = coeffs.subs({dt: di, it: shifts[i][di]})
            # Track Indexed substitutions
            idx_subs.append(mapper)

        # Temporaries for the indirection dimensions
        temps = [Eq(v, k, implicit_dims=implicit_dims) for k, v in points.items()]

        return idx_subs, temps

    @property
    def _interpolation_coeffs(self):
        return self.sfunction.interpolation_coeffs

    @property
    def _interpolation_coeffsp(self):
        d = self.sfunction.interpolation_coeffs.dimensions[1]
        return prod([self.sfunction.interpolation_coeffs._subs(d, i)
                     for (i, _) in enumerate(self.sfunction.grid.dimensions)])

    def subs_coords(self, _expr, *idx_subs):
        b = self._interpolation_coeffsp
        return [_expr.xreplace(v_sub) * b.xreplace(v_sub) for v_sub in idx_subs]

    def subs_coords_eq(self, field, _expr, *idx_subs, implicit_dims=None):
        b = self._interpolation_coeffsp
        return [Inc(field.xreplace(vsub), _expr.xreplace(vsub) * b.xreplace(vsub),
                    implicit_dims=implicit_dims) for vsub in idx_subs]
