from devito import *
from examples.seismic.source import WaveletSource, RickerSource, GaborSource, TimeAxis
from examples.seismic import plot_image
import numpy as np
from devito.tools import as_tuple

shape = (101, 101, 101)
nt = 30
so = 4

# Initial grid: km x km, with spacing
shape = shape  # Number of grid point (nx, nz)
spacing = as_tuple(10.0 for _ in range(len(shape)))
extent = tuple([s*(n-1) for s, n in zip(spacing, shape)])

x = SpaceDimension(name='x', spacing=Constant(name='h_x', value=extent[0]/(shape[0]-1)))  # noqa
y = SpaceDimension(name='y', spacing=Constant(name='h_y', value=extent[1]/(shape[1]-1)))  # noqa
z = SpaceDimension(name='z', spacing=Constant(name='h_z', value=extent[2]/(shape[2]-1)))  # noqa
grid = Grid(extent=extent, shape=shape, dimensions=(x, y, z))

# Timestep size from Eq. 7 with V_p=6000. and dx=100
t0, tn = 0., nt
dt = (10. / np.sqrt(2.)) / 6.
time_range = TimeAxis(start=t0, stop=tn, step=dt)

src = RickerSource(name='src', grid=grid, f0=0.01, time_range=time_range)
src.coordinates.data[:] = np.array([250., 250., 250.])

# Now we create the velocity and pressure fields
v = VectorTimeFunction(name='v', grid=grid, space_order=so, time_order=1)
tau = TensorTimeFunction(name='t', grid=grid, space_order=so, time_order=1)

# We need some initial conditions
V_p = 2.0
V_s = 1.0
density = 1.8

# The source injection term
src_xx = src.inject(field=tau.forward[0, 0], expr=src)
src_yy = src.inject(field=tau.forward[1, 1], expr=src)
src_zz = src.inject(field=tau.forward[2, 2], expr=src)

# Thorbecke's parameter notation
cp2 = V_p*V_p
cs2 = V_s*V_s
ro = 1/density

mu = cs2*density
l = (cp2*density - 2*mu)

# First order elastic wave equation
pde_v = v.dt - ro * div(tau)
pde_tau = (tau.dt - l * diag(div(v.forward)) - mu * (grad(v.forward) +
           grad(v.forward).transpose(inner=False)))

# Time update
u_v = Eq(v.forward, solve(pde_v, v.forward))
u_t = Eq(tau.forward, solve(pde_tau, tau.forward))

# Inject sources. We use it to preinject data
op = Operator([u_v] + [u_t] + src_xx + src_yy + src_zz)
op(dt=dt)

dv_norm_v0 = norm(v[0])
dv_norm_v1 = norm(v[1])
dv_norm_tau0 = norm(tau[0])
dv_norm_tau1 = norm(tau[1])
dv_norm_tau2 = norm(tau[2])
dv_norm_tau3 = norm(tau[3])

print(dv_norm_v0)
print(dv_norm_v1)
print(dv_norm_tau0)
print(dv_norm_tau1)
print(dv_norm_tau2)
print(dv_norm_tau3)
