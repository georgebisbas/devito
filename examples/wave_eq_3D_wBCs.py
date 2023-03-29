# A 3D wave equation using Devito
# BC modelling included
# PyVista plotting included

import argparse
import numpy as np

from devito import Grid, TimeFunction, Eq, solve, Operator, Constant, norm
from examples.seismic import Model

parser = argparse.ArgumentParser(description='Process arguments.')

parser.add_argument("-d", "--shape", default=(11, 11, 11), type=int, nargs="+",
                    help="Number of grid points along each axis")
parser.add_argument("-so", "--space_order", default=4,
                    type=int, help="Space order of the simulation")
parser.add_argument("-to", "--time_order", default=1,
                    type=int, help="Time order of the simulation")
parser.add_argument("-nt", "--nt", default=40,
                    type=int, help="Simulation time in millisecond")
parser.add_argument("-bls", "--blevels", default=2, type=int, nargs="+",
                    help="Block levels")
parser.add_argument("-plot", "--plot", default=False, type=bool, help="Plot3D")
args = parser.parse_args()


def plot_3dfunc(u):
    # Plot a 3D structured grid using pyvista

    import matplotlib.pyplot as plt
    import pyvista as pv

    cmap = plt.colormaps["viridis"]
    values = u.data[0, :, :, :]
    vistagrid = pv.UniformGrid()
    vistagrid.dimensions = np.array(values.shape) + 1
    vistagrid.spacing = (1, 1, 1)
    vistagrid.origin = (0, 0, 0)  # The bottom left corner of the data set
    vistagrid.cell_data["values"] = values.flatten(order="F")
    vistaslices = vistagrid.slice_orthogonal()
    # vistagrid.plot(show_edges=True)
    vistaslices.plot(cmap=cmap)


# Some variable declarations
nx, ny, nz = args.shape
nt = args.nt
nu = .5
dx = 2. / (nx - 1)
dy = 2. / (ny - 1)
dz = 2. / (nz - 1)
sigma = .25
dt = sigma * dx * dz * dy / nu
so = args.space_order
to = 2
print("dt %s, dx %s, dy %s, dz %s" % (dt, dx, dy, dz))


# Define a physical size
shape = args.shape  # Number of grid point (nx, nz)
spacing = (10., 10., 10.)  # Grid spacing in m. The domain size is now 1km by 1km
origin = (0., 0., 0.)  # What is the location of the top left corner. This is necessary to define
# the absolute location of the source and receivers

grid = Grid(shape=(nx, ny, nz), extent=(2., 2., 2.))
u = TimeFunction(name='u', grid=grid, space_order=so, time_order=2)
# init_hat(field=u.data[0], dx=dx, dy=dy, value=2.)
# u.data[:, :, :, :] = 1.5
# u.data[:, int(nx/2), int(ny/2), int(nz/2)] = 5

# Create an equation with second-order derivatives

# We can now write the PDE
pde = 0.1*u.dt2 - u.laplace + 0.2*u.dt

# The PDE representation is as on paper
pde

stencil = Eq(u.forward, solve(pde, u.forward))
stencil


# Create boundary condition expressions
x, y, z = grid.dimensions
t = grid.stepping_dim

# Add boundary conditions
bc = [Eq(u[t+1, x, y, 0], 2.)]  # bottom
bc += [Eq(u[t+1, x, y, nz-1], 2.)]  # top

bc += [Eq(u[t+1, 0, y, z], 2.)]  # left
bc += [Eq(u[t+1, nx-1, y, z], 2.)]  # right

bc += [Eq(u[t+1, x, 0, z], 2.)]  # front
bc += [Eq(u[t+1, x, ny-1, z], 2.)]  # back

print(stencil)

# Create an operator that updates the forward stencil point
# plus adding boundary conditions
op = Operator([stencil] + bc, subdomain=grid.interior)

# No BCs
# op = Operator([stencil], subdomain=grid.interior)
print(op.ccode)

# Apply the operator for a number of timesteps
op(time=nt, dt=dt)

if args.plot:
    plot_3dfunc(u)

print("Field norm is:", norm(u))
import pdb;pdb.set_trace()