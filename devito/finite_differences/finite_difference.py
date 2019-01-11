from functools import partial, wraps

from sympy import S, finite_diff_weights

from devito.finite_differences import Differentiable
from devito.tools import Tag

__all__ = ['first_derivative', 'second_derivative', 'cross_derivative',
           'generic_derivative', 'left', 'right', 'centered', 'transpose']

# Number of digits for FD coefficients to avoid roundup errors and non-deterministic
# code generation
_PRECISION = 9


class Transpose(Tag):
    """
    Utility class to change the sign of a derivative. This is only needed
    for odd order derivatives, which require a minus sign for the transpose.
    """
    pass


direct = Transpose('direct', 1)
transpose = Transpose('transpose', -1)


class Side(Tag):
    """
    Class encapsulating the side of the shift for derivatives.
    """

    def adjoint(self, matvec):
        if matvec == direct:
            return self
        else:
            if self == centered:
                return centered
            elif self == right:
                return left
            elif self == left:
                return right
            else:
                raise ValueError("Unsupported side value")


left = Side('left', -1)
right = Side('right', 1)
centered = Side('centered', 0)


def check_input(func):
    @wraps(func)
    def wrapper(expr, *args, **kwargs):
        if expr.is_Number:
            return S.Zero
        elif not isinstance(expr, Differentiable):
            raise ValueError("`%s` must be of type Differentiable (found `%s`)"
                             % (expr, type(expr)))
        else:
            return func(expr, *args, **kwargs)
    return wrapper


@check_input
def first_derivative(expr, **kwargs):
    """
    First-order derivative of a given expression.

    Parameters
    ----------
    expr : expr-like
        Expression for which the first-order derivative is produced.
    **kwargs
        - ``dim``: Dimension w.r.t. which to differentiate.
        - ``diff``: Finite-difference symbol, defaults to `h`.
        - ``order``: Discretisation order of the derivative.

    Returns
    -------
    expr-like
        The first-order derivative of ``expr``

    Examples
    --------
    >>> from devito import Function, Grid, first_derivative
    >>> grid = Grid(shape=(4, 4))
    >>> x, _ = grid.dimensions
    >>> f = Function(name='f', grid=grid)
    >>> g = Function(name='g', grid=grid)
    >>> first_derivative(f*g, dim=x)
    -f(x, y)*g(x, y)/h_x + f(x + h_x, y)*g(x + h_x, y)/h_x

    This is also more easily obtainable via:

    >>> (f*g).dx
    -f(x, y)*g(x, y)/h_x + f(x + h_x, y)*g(x + h_x, y)/h_x
    """
    dim = kwargs.get('dim')
    diff = kwargs.get('diff', dim.spacing)
    order = int(kwargs.get('order', 1))
    matvec = kwargs.get('matvec', direct)
    side = kwargs.get('side', centered).adjoint(matvec)
    deriv = 0
    # Stencil positions for non-symmetric cross-derivatives with symmetric averaging
    if side == right:
        ind = [(dim + i * diff) for i in range(-int(order / 2) + 1 - (order % 2),
                                               int((order + 1) / 2) + 2 - (order % 2))]
    elif side == left:
        ind = [(dim - i * diff) for i in range(-int(order / 2) + 1 - (order % 2),
                                               int((order + 1) / 2) + 2 - (order % 2))]
    else:
        ind = [(dim + i * diff) for i in range(-int(order / 2),
                                               int((order + 1) / 2) + 1)]
    # Finite difference weights from Taylor approximation with this positions
    c = finite_diff_weights(1, ind, dim)
    c = c[-1][-1]
    all_dims = tuple(set((dim,) + tuple([i for i in expr.indices if i.root == dim])))
    # Loop through positions
    for i in range(0, len(ind)):
            subs = dict([(d, ind[i].subs({dim: d})) for d in all_dims])
            deriv += expr.subs(subs) * c[i]
    return (matvec.val*deriv).evalf(_PRECISION)


@check_input
def second_derivative(expr, dim, fd_order, stagger=None):
    """
    Second-order derivative of a given expression.

    Parameters
    ----------
    expr : expr-like
        Expression for which the derivative is produced.
    dim : Dimension
        The Dimension w.r.t. which to differentiate.
    fd_order : int
        Coefficient discretization order. Note: this impacts the width of
        the resulting stencil.
    stagger : Side, optional
        Shift of the finite-difference approximation.

    Returns
    -------
    expr-like
        The ``deriv-order`` orderderivative of ``expr``.

    Examples
    --------
    >>> from devito import Function, Grid, second_derivative
    >>> grid = Grid(shape=(4, 4))
    >>> x, _ = grid.dimensions
    >>> f = Function(name='f', grid=grid, space_order=2)
    >>> g = Function(name='g', grid=grid, space_order=2)
    >>> second_derivative(f*g, dim=x, fd_order=2)
    -2.0*f(x, y)*g(x, y)/h_x**2 + f(x - h_x, y)*g(x - h_x, y)/h_x**2 +\
 f(x + h_x, y)*g(x + h_x, y)/h_x**2

    This is also more easily obtainable via:

    >>> (f*g).dx2
    -2.0*f(x, y)*g(x, y)/h_x**2 + f(x - h_x, y)*g(x - h_x, y)/h_x**2 +\
 f(x + h_x, y)*g(x + h_x, y)/h_x**2
    """

    return generic_derivative(expr, 2, dim, fd_order, stagger=None)


@check_input
def cross_derivative(expr, dims, fd_order, deriv_order, stagger=(None, None)):
    """
    Arbitrary-order cross derivative of a given expression.

    Parameters
    ----------
    expr : expr-like
        Expression for which the cross derivative is produced.
    dims : tuple
        Dimensions w.r.t. which to differentiate.
    fd_order : int
        Coefficient discretization order. Note: this impacts the width of
        the resulting stencil.
    deriv_order : int
        Derivative order, e.g. 2 for a second-order derivative.

    Examples
    --------
    >>> from devito import Function, Grid, second_derivative
    >>> grid = Grid(shape=(4, 4))
    >>> x, y = grid.dimensions
    >>> f = Function(name='f', grid=grid, space_order=2)
    >>> g = Function(name='g', grid=grid, space_order=2)
    >>> cross_derivative(f*g, dims=(x, y), fd_order=(2, 2), deriv_order=(1, 1))
    -0.5*(-0.5*f(x - h_x, y - h_y)*g(x - h_x, y - h_y)/h_x +\
 0.5*f(x + h_x, y - h_y)*g(x + h_x, y - h_y)/h_x)/h_y +\
 0.5*(-0.5*f(x - h_x, y + h_y)*g(x - h_x, y + h_y)/h_x +\
 0.5*f(x + h_x, y + h_y)*g(x + h_x, y + h_y)/h_x)/h_y

    This is also more easily obtainable via:

    >>> (f*g).dxdy
    -0.5*(-0.5*f(x - h_x, y - h_y)*g(x - h_x, y - h_y)/h_x +\
 0.5*f(x + h_x, y - h_y)*g(x + h_x, y - h_y)/h_x)/h_y +\
 0.5*(-0.5*f(x - h_x, y + h_y)*g(x - h_x, y + h_y)/h_x +\
 0.5*f(x + h_x, y + h_y)*g(x + h_x, y + h_y)/h_x)/h_y
    """
    first = expr
    for d, fd, dim, s in zip(deriv_order, fd_order, dims, stagger):
        first = generic_derivative(first, dim=dim, fd_order=fd, deriv_order=d, stagger=s)

    return first


@check_input
def generic_derivative(expr, dim, fd_order, deriv_order, stagger=None):
    """
    Arbitrary-order derivative of a given expression.

    Parameters
    ----------
    expr : expr-like
        Expression for which the derivative is produced.
    deriv_order : int
        Derivative order, e.g. 2 for a second-order derivative.
    dim : Dimension
        The Dimension w.r.t. which to differentiate.
    fd_order : int
        Coefficient discretization order. Note: this impacts the width of
        the resulting stencil.
    stagger : Side, optional
        Shift of the finite-difference approximation.

    Returns
    -------
    expr-like
        The derivative of ``expr`` of order ``deriv-order``.
    """

    diff = dim.spacing

    if stagger == left or not expr.is_Staggered:
        off = -.5
    elif stagger == right:
        off = .5
    else:
        off = 0

    if expr.is_Staggered:
        indices = list(set([(dim + int(i+.5+off) * dim.spacing)
                            for i in range(-fd_order//2, fd_order//2)]))
        x0 = (dim + off*diff)
        if fd_order < 2:
            indices = [dim + diff, dim] if stagger == right else [dim - diff, dim]

    else:
        indices = [(dim + i * dim.spacing) for i in range(-fd_order//2, fd_order//2 + 1)]
        x0 = dim
        if fd_order < 2:
            indices = [dim, dim + diff]

    c = finite_diff_weights(deriv_order, indices, x0)[-1][-1]

    deriv = 0
    all_dims = tuple(set((dim, ) +
                     tuple([i for i in expr.indices if i.root == dim])))
    for i in range(0, len(indices)):
            subs = dict([(d, indices[i].subs({dim: d})) for d in all_dims])
            deriv += expr.subs(subs) * c[i]

    return deriv.evalf(_PRECISION)


def generate_fd_shortcuts(function):
    """Create all legal finite-difference derivatives for the given Function."""
    dimensions = function.indices
    space_fd_order = function.space_order
    time_fd_order = function.time_order if (function.is_TimeFunction or
                                            function.is_SparseTimeFunction) else 0

    deriv_function = generic_derivative
    c_deriv_function = cross_derivative

    side = dict()
    for (d, s) in zip(dimensions, function.staggered):
        if s == 0:
            side[d] = left
        elif s == 1:
            side[d] = right
        else:
            side[d] = centered

    derivatives = dict()
    done = []
    for d in dimensions:
        # Dimension is treated, remove from list
        done += [d]
        other_dims = tuple(i for i in dimensions if i not in done)
        # Dimension name and corresponding FD order
        dim_order = time_fd_order if d.is_Time else space_fd_order
        name = 't' if d.is_Time else d.root.name
        # All possible derivatives go up to the dimension FD order
        for o in range(1, dim_order + 1):
            deriv = partial(deriv_function, deriv_order=o, dim=d,
                            fd_order=dim_order, stagger=side[d])
            name_fd = 'd%s%d' % (name, o) if o > 1 else 'd%s' % name
            desciption = 'derivative of order %d w.r.t dimension %s' % (o, d)

            derivatives[name_fd] = (deriv, desciption)
            # Cross derivatives with the other dimension
            # Skip already done dimensions a dxdy is the same as dydx
            for d2 in other_dims:
                dim_order2 = time_fd_order if d2.is_Time else space_fd_order
                name2 = 't' if d2.is_Time else d2.root.name
                for o2 in range(1, dim_order2 + 1):
                    deriv = partial(c_deriv_function, deriv_order=(o, o2), dims=(d, d2),
                                    fd_order=(dim_order, dim_order2),
                                    stagger=(side[d], side[d2]))
                    name_fd2 = 'd%s%d' % (name, o) if o > 1 else 'd%s' % name
                    name_fd2 += 'd%s%d' % (name2, o2) if o2 > 1 else 'd%s' % name2
                    desciption = 'derivative of order (%d, %d) ' % (o, o2)
                    desciption += 'w.r.t dimension (%s, %s) ' % (d, d2)
                    derivatives[name_fd2] = (deriv, desciption)

    # Add non-conventional, non-centered first-order FDs
    for d in dimensions:
        name = 't' if d.is_Time else d.root.name
        if function.is_Staggered:
            # Add centered first derivatives if staggered
            deriv = partial(deriv_function, deriv_order=1, dim=d,
                            fd_order=dim_order, stagger=centered)
            name_fd = 'd%sc' % name
            desciption = 'centered derivative staggered w.r.t dimension %s' % d

            derivatives[name_fd] = (deriv, desciption)
        else:
            # Left
            dim_order = time_fd_order if d.is_Time else space_fd_order
            deriv = partial(first_derivative, order=dim_order, dim=d, side=left)
            name_fd = 'd%sl' % name
            desciption = 'left first order derivative w.r.t dimension %s' % d
            derivatives[name_fd] = (deriv, desciption)
            # Right
            deriv = partial(first_derivative, order=dim_order, dim=d, side=right)
            name_fd = 'd%sr' % name
            desciption = 'right first order derivative w.r.t dimension %s' % d
            derivatives[name_fd] = (deriv, desciption)

    return derivatives
