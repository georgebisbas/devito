from devito import *
from examples.seismic.source import (RickerSource, TimeAxis, Receiver)
from examples.seismic import plot_image
import numpy as np

from sympy import init_printing, latex
init_printing(use_latex='mathjax')

# Initial grid: 1km x 1km, with spacing 100m
extent = (1500., 1500.)
shape = (201, 201)
x = SpaceDimension(name='x', spacing=Constant(name='h_x', value=extent[0]/(shape[0]-1)))
z = SpaceDimension(name='z', spacing=Constant(name='h_z', value=extent[1]/(shape[1]-1)))
grid = Grid(extent=extent, shape=shape, dimensions=(x, z))


# Timestep size from Eq. 7 with V_p=6000. and dx=100
t0, tn = 0., 300.
dt = (10. / np.sqrt(2.)) / 6.
time_range = TimeAxis(start=t0, stop=tn, step=dt)

src = RickerSource(name='src', grid=grid, f0=0.01, time_range=time_range)

src.coordinates.data[:] = np.array([250., 250.])
# Now we create the velocity and pressure fields
so = 2

v = VectorTimeFunction(name='v', grid=grid, space_order=so, time_order=1)
tau = TensorTimeFunction(name='t', grid=grid, space_order=so, time_order=1)

# Now let's try and create the staggered updates
t = grid.stepping_dim
time = grid.time_dim

# The source injection term
src_xx = src.inject(field=tau.forward[0, 0], expr=src)

# The receiver
nrec = 1
rec = Receiver(name="rec", grid=grid, npoint=nrec, time_range=time_range)
rec.coordinates.data[:, 0] = np.linspace(0., extent[0], num=nrec)
rec.coordinates.data[:, -1] = 5.

rec_term = rec.interpolate(expr=v[0])

# First order elastic wave equation
pde_v = v.dt - div(tau)
pde_tau = (tau.dt - (grad(v.forward)))

# Time update
u_v = Eq(v.forward, solve(pde_v, v.forward))
u_t = Eq(tau.forward, solve(pde_tau, tau.forward))

op = Operator([u_v] + [u_t] + src_xx + rec_term)

op(dt=dt)
print(op.ccode)

# import pdb; pdb.set_trace();

print(norm(v[0]))
plot_image(v[0].data[0])
