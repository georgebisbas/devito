import numpy as np

from examples.cfd import plot_field, init_smooth, init_hat
from devito import Grid, TimeFunction, Eq, solve, Operator, Constant, norm
from matplotlib.pyplot import pause

from sympy import nsimplify, pprint

# Some variable declarations
nx = 80
ny = 80
nt = 500
nu = .5
dx = 2. / (nx - 1)
dy = 2. / (ny - 1)
sigma = .25
dt = sigma * dx * dy / nu
print("dx %s, dy %s" % (dx, dy))

grid = Grid(shape=(nx, ny), extent=(2., 2.))
u = TimeFunction(name='u', grid=grid, space_order=2)
# init_hat(field=u.data[0], dx=dx, dy=dy, value=2.)
u.data[:, :, :] = 1
init_smooth(field=u.data[0], dx=dx, dy=dy)

# Note u.data[0] == u.data[0,:,:]

a = Constant(name='a')
# Create an equation with second-order derivatives
eq = Eq(u.dt, a * u.laplace, subdomain=grid.interior)
stencil = solve(eq, u.forward)
eq_stencil = Eq(u.forward, stencil)

# Create boundary condition expressions
x, y = grid.dimensions
t = grid.stepping_dim
bc = [Eq(u[t+1, 0, y], 1.)]  # left
bc += [Eq(u[t+1, nx-1, y], 2.)]  # right
bc += [Eq(u[t+1, x, ny-1], 1.)]  # top
bc += [Eq(u[t+1, x, 0], 1.)]  # bottom

print(eq_stencil)

plot_field(u.data[0], zmax=4); pause(1)

# Create an operator that updates the forward stencil point
op = Operator([eq_stencil] + bc, subdomain=grid.interior)

# Apply the operator for a number of timesteps
op(time=nt, dt=dt, a=nu)

# nsimplify: 1.0*x = x
# pprint(nsimplify(stencil))

plot_field(u.data[0], zmax=4); pause(1)


import pdb;pdb.set_trace()
