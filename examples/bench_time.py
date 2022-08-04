
import pytest
import numpy as np

from devito import (Grid, Eq, Function, TimeFunction, Operator, norm,
                    Constant, solve)
from devito.ir import Expression, Iteration, FindNodes, FindSymbols
import argparse

parser = argparse.ArgumentParser(description='Process arguments.')

parser.add_argument("-d", "--shape", default=(11, 11, 11), type=int, nargs="+",
                    help="Number of grid points along each axis")
parser.add_argument("-so", "--space_order", default=4,
                    type=int, help="Space order of the simulation")
parser.add_argument("-to", "--time_order", default=2,
                    type=int, help="Time order of the simulation")
parser.add_argument("-tn", "--tn", default=64,
                    type=int, help="Simulation time in millisecond")
parser.add_argument("-bls", "--blevels", default=2, type=int, nargs="+",
                    help="Block levels")
args = parser.parse_args()


nx, ny, nz = args.shape
tn = args.tn
nu = .5
dx = 2. / (nx - 1)
dy = 2. / (ny - 1)
dz = 2. / (nz - 1)
sigma = .25
dt = sigma * dx * dz * dy / nu


so = args.space_order
to = 1

# Initialise u
init_value = 6.5

# Field initialization
grid = Grid(shape=(nx, ny, nz))
u = TimeFunction(name='u', grid=grid, space_order=so, time_order=to)
u.data[:, :, :] = init_value

# Create an equation with second-order derivatives
a = Constant(name='a')
eq = Eq(u.dt, a*u.laplace + 0.1, subdomain=grid.interior)
stencil = solve(eq, u.forward)
eq0 = Eq(u.forward, stencil)

# List comprehension would need explicit locals/globals mappings to eval
op0 = Operator(eq0, opt=('advanced'))
# op0.apply(time_M=tn, dt=dt, **{'x0_blk0_size': 32, 'y0_blk0_size': 16})
# op0.apply(time_M=tn, dt=dt, **{'x0_blk0_size': 64, 'y0_blk0_size': 32})
op0.apply(time_M=tn-1, dt=dt)
norm_u = norm(u)
print(norm_u)


u.data[:] = init_value
#
op1 = Operator(eq0, opt=('advanced', {'skewing': True,
                         'blocklevels': args.blevels}))
# op1.apply(time_M=tn, dt=dt, **{'time0_blk0_size': 64, 'x0_blk0_size': 32, 'x0_blk1_size': 8, 'y0_blk0_size': 32, 'y0_blk1_size': 4})   # Large problem 1024^3 - 256

# op1.apply(time_M=tn, dt=dt, **{'time0_blk0_size': 64, 'x0_blk0_size': 32, 'x0_blk1_size': 4, 'y0_blk0_size': 32, 'y0_blk1_size': 4})   # Large problem 1024^3 - 256
# op1.apply(time_M=tn, dt=dt, **{'time0_blk0_size': 64, 'x0_blk0_size': 64, 'x0_blk1_size': 16, 'y0_blk0_size': 32, 'y0_blk1_size': 8})  # Medium problem 512^3 - 256
# op1.apply(time_M=tn, dt=dt, **{'time0_blk0_size': 64, 'x0_blk0_size': 64, 'x0_blk1_size': 16, 'y0_blk0_size': 32, 'y0_blk1_size': 8})  # Medium problem 512^3 - 256

# op1.apply(time_M=tn, dt=dt, **{'time0_blk0_size': 64, 'x0_blk0_size': 64, 'x0_blk1_size': 4, 'y0_blk0_size': 32, 'y0_blk1_size': 4})  # Medium problem 512^3 - 512 gcc

#op1.apply(time_M=tn, dt=dt, **{'time0_blk0_size': 64, 'x0_blk0_size': 64, 'x0_blk1_size': 4, 'y0_blk0_size': 64, 'y0_blk1_size': 4})  # Medium problem 512^3 - 512 gcc

# op1.apply(time_M=tn, dt=dt, **{'time0_blk0_size': 64, 'x0_blk0_size': 32, 'x0_blk1_size': 8, 'y0_blk0_size': 32, 'y0_blk1_size': 8})  # Medium problem 512^3 - 256
# op1.apply(time_M=tn, dt=dt, **{'time0_blk0_size': 512, 'x0_blk0_size': 2112, 'x0_blk1_size': 24, 'y0_blk0_size': 2112, 'y0_blk1_size': 24})  # Medium problem 512^3 - 256
op1.apply(time_M=tn, dt=dt)

print(norm_u)

print(norm(u))
# assert np.isclose(norm(u), norm_u, atol=1e-4, rtol=0)

# iters = FindNodes(Iteration).visit(op1)