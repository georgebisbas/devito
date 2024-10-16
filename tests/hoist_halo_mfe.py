from devito import *
from devito.ir.iet import FindNodes, Call

grid = Grid(shape=(12,))
x = grid.dimensions[0]
t = grid.stepping_dim

i = Dimension(name='i')

f = TimeFunction(name='f', grid=grid)
g = Function(name='g', grid=grid)
h = Function(name='h', grid=grid)

op = Operator([Eq(f.forward, f[t, x-1] + f[t, x+1] + 1.),
               Inc(g[i], f[t, h[i]] + 1.)])

calls = FindNodes(Call).visit(op)
assert len(calls) == 1

# Below, there is a flow-dependence along `x`, so a further halo update
# before the Inc is required
op = Operator([Eq(f.forward, f[t, x-1] + f[t, x+1] + 1.),
               Inc(g[i], f[t+1, h[i]] + 1.)])

op.apply(i_M=4, time_M=2)

calls = FindNodes(Call).visit(op)
import pdb;pdb.set_trace()
assert len(calls) == 2


