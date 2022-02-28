from sympy import sympify
from collections import Counter
from sympy import Mod

from devito.ir.clusters import Queue
from devito.ir.support import (AFFINE, PARALLEL, PARALLEL_IF_ATOMIC, PARALLEL_IF_PVT,
                               SEQUENTIAL, SKEWABLE, TILABLE, Interval, IntervalGroup,
                               IterationSpace, Scope)
from devito.symbolics import uxreplace, INT, xreplace_indices, evalrel, retrieve_indexed
from devito.symbolics import uxreplace, xreplace_indices
from devito.tools import UnboundedMultiTuple, as_tuple, flatten
from devito.types import BlockDimension


from devito.types import RIncrDimension
from devito.tools import as_list, as_tuple

__all__ = ['blocking', 'skewing']


def blocking(clusters, sregistry, options):
    """
    Loop blocking to improve data locality.

    Parameters
    ----------
    clusters : tuple of Clusters
        Input Clusters, subject of the optimization pass.
    options : dict
        The optimization options.
        * `blockinner` (boolean, False): enable/disable loop blocking along the
           innermost loop.
        * `blocklevels` (int, 1): 1 => classic loop blocking; 2 for two-level
           hierarchical blocking.

    Example
    -------
        * A typical use case, e.g.
          .. code-block::
                            Classical   +blockinner  2-level Hierarchical
            for x            for xb        for xb         for xbb
              for y    -->    for yb        for yb         for ybb
                for z          for x         for zb         for xb
                                for y         for x          for yb
                                 for z         for y          for x
                                                for z          for y
                                                                for z
    """
    if options['blockrelax']:
        analyzer = AnalyzeBlocking()
    else:
        analyzer = AnalyzeHeuristicBlocking(options)
    clusters = analyzer.process(clusters)
    clusters = AnalyzeSkewing().process(clusters)

    if options['blocklevels'] > 0:
        clusters = SynthesizeBlocking(sregistry, options).process(clusters)

    return clusters


class AnayzeBlockingBase(Queue):

    """
    Encode the TILABLE property.
    """

    def process(self, clusters):
        return self._process_fatd(clusters, 1)

    def _process_fatd(self, clusters, level, prefix=None):
        # Truncate recursion in case of TILABLE, non-perfect sub-nests, as
        # it's an unsupported case
        if prefix:
            d = prefix[-1].dim

            if any(TILABLE in c.properties[d] for c in clusters) and \
               len({c.itintervals[:level] for c in clusters}) > 1:
                return clusters

        return super()._process_fatd(clusters, level, prefix)


class AnalyzeBlocking(AnayzeBlockingBase):

    def callback(self, clusters, prefix):
        if not prefix:
            return clusters

        d = prefix[-1].dim

        for c in clusters:
            if not {PARALLEL,
                    PARALLEL_IF_ATOMIC,
                    PARALLEL_IF_PVT}.intersection(c.properties[d]):
                return clusters

        # All good, `d` is actually TILABLE
        processed = attach_property(clusters, d, TILABLE)

        return processed


class AnalyzeHeuristicBlocking(AnayzeBlockingBase):

    def __init__(self, options):
        super().__init__()

        self.inner = options['blockinner']

    def process(self, clusters):
        clusters = super().process(clusters)

        # Heuristic: if there aren't at least two TILABLE Dimensions, drop it
        processed = []
        for c in clusters:
            ntilable = len([TILABLE for v in c.properties.values() if TILABLE in v])
            if ntilable > 1:
                processed.append(c)
            else:
                properties = {d: v - {TILABLE} for d, v in c.properties.items()}
                processed.append(c.rebuild(properties=properties))

        return processed

    def callback(self, clusters, prefix):
        if not prefix:
            return clusters

        d = prefix[-1].dim

        for c in clusters:
            # PARALLEL* and AFFINE are necessary conditions
            if AFFINE not in c.properties[d] or \
               not ({PARALLEL, PARALLEL_IF_PVT} & c.properties[d]):
                return clusters

            # Heuristic: innermost Dimensions may be ruled out a-priori
            is_inner = d is c.itintervals[-1].dim
            if is_inner and not self.inner:
                return clusters

            # Heuristic: TILABLE not worth it if not within a SEQUENTIAL Dimension
            if not any(SEQUENTIAL in c.properties[i.dim] for i in prefix[:-1]):
                return clusters

            # Heuristic: same as above if there's a local SubDimension
            if any(i.dim.is_Sub and i.dim.local for i in c.itintervals):
                return clusters

        if len(clusters) > 1:
            # Heuristic: same as above if it induces dynamic bounds
            exprs = flatten(c.exprs for c in as_tuple(clusters))
            scope = Scope(exprs)
            if any(i.is_lex_non_stmt for i in scope.d_all_gen()):
                return clusters
        else:
            # Just avoiding potentially expensive checks
            pass

        # All good, `d` is actually TILABLE
        processed = attach_property(clusters, d, TILABLE)

        return processed


class AnalyzeSkewing(Queue):

    """
    Encode the SKEWABLE Dimensions.
    """

    def callback(self, clusters, prefix):
        if not prefix:
            return clusters

        d = prefix[-1].dim

        for c in clusters:
            if TILABLE not in c.properties[d]:
                return clusters

        processed = attach_property(clusters, d, SKEWABLE)

        return processed


class SynthesizeBlocking(Queue):

    template = "%s%d_blk%s"

    def __init__(self, sregistry, options):
        self.sregistry = sregistry

        self.levels = options['blocklevels']

        # A tool to unroll the explicit integer block shapes, should there be any
        if options['par-tile']:
            self.blk_size_gen = UnboundedMultiTuple(*options['par-tile'])
        else:
            self.blk_size_gen = None

        super().__init__()

    def process(self, clusters):
        return self._process_fatd(clusters, 1)

    def _make_key_hook(self, cluster, level):
        return (tuple(cluster.guards.get(i.dim) for i in cluster.itintervals[:level]),)

    def callback(self, clusters, prefix):
        if not prefix:
            return clusters

        d = prefix[-1].dim

        if not any(TILABLE in c.properties[d] for c in clusters):
            return clusters

        # Create the block Dimensions (in total `self.levels` Dimensions)
        base = self.sregistry.make_name(prefix=d.name)

        if self.blk_size_gen:
            # If a new TILABLE nest, pull what would be the next par-tile entry
            if not any(i.dim.is_Block for i in prefix):
                self.blk_size_gen.iter()

            step = sympify(self.blk_size_gen.next())
        else:
            # This will result in a parametric step, e.g. `x0_blk0_size`
            step = None

        name = self.sregistry.make_name(prefix="%s_blk" % base)
        bd = BlockDimension(name, d, d.symbolic_min, d.symbolic_max, step)
        step = bd.step
        block_dims = [bd]

        for _ in range(1, self.levels):
            name = self.sregistry.make_name(prefix="%s_blk" % base)
            bd = BlockDimension(name, bd, bd, bd + bd.step - 1, size=step)
            block_dims.append(bd)

        bd = BlockDimension(d.name, bd, bd, bd + bd.step - 1, 1, size=step)
        block_dims.append(bd)

        # name = self.template % (d.name, self.nblocked[d], '%d')
        # block_dims = create_block_dims(name, d, self.levels)

        processed = []
        for c in clusters:
            if TILABLE in c.properties[d]:
                ispace = decompose(c.ispace, d, block_dims)

                # Use the innermost IncrDimension in place of `d`
                exprs = [uxreplace(e, {d: block_dims[-1]}) for e in c.exprs]

                # The new Cluster properties
                # TILABLE property is dropped after the blocking.
                properties = dict(c.properties)
                properties.pop(d)
                properties.update({bd: c.properties[d] - {TILABLE} for bd in block_dims})

                processed.append(c.rebuild(exprs=exprs, ispace=ispace,
                                           properties=properties))
            else:
                processed.append(c)

        return processed


def preprocess(clusters, options):
    # Preprocess: heuristic: drop TILABLE from innermost Dimensions to
    # maximize vectorization
    inner = bool(options['blockinner'])
    processed = []
    for c in clusters:
        ntilable = len([i for i in c.properties.values() if TILABLE in i])
        ntilable -= int(not inner)
        if ntilable <= 1:
            properties = {k: v - {TILABLE, SKEWABLE} for k, v in c.properties.items()}
            processed.append(c.rebuild(properties=properties))
        elif not inner:
            d = c.itintervals[-1].dim
            properties = dict(c.properties)
            properties[d] = properties[d] - {TILABLE, SKEWABLE}
            processed.append(c.rebuild(properties=properties))
        else:
            processed.append(c)

    return processed


def create_block_dims(name, d, levels, **kwargs):
    """
    Create the block Dimensions (in total `self.levels` Dimensions)
    """
    sf = kwargs.pop('sf', 1)
    bd = RIncrDimension(name % 0, d, d.symbolic_min, d.symbolic_max,
                        rmax=sf*d.symbolic_max)
    size = bd.step
    block_dims = [bd]

    for i in range(1, levels):
        bd = RIncrDimension(name % i, bd, bd, bd + bd.step - 1, size=size)
        block_dims.append(bd)

    bd = RIncrDimension(d.name, bd, bd, bd + bd.step - 1, 1, size=size,
                        rmax=evalrel(min, [bd + bd.step - 1, sf*d.root.symbolic_max]),
                        rstep=sf)
    block_dims.append(bd)

    return block_dims


def decompose(ispace, d, block_dims):
    """
    Create a new IterationSpace in which the `d` Interval is decomposed
    into a hierarchy of Intervals over ``block_dims``.
    """
    # Create the new Intervals
    intervals = []
    for i in ispace:
        if i.dim is d:
            intervals.append(i.switch(block_dims[0]))
            intervals.extend([i.switch(bd).zero() for bd in block_dims[1:]])
        else:
            intervals.append(i)

    # Create the intervals relations
    # 1: `bd > d`
    relations = [tuple(block_dims)]

    # 2: Suitably replace `d` with all `bd`'s
    for r in ispace.relations:
        try:
            n = r.index(d)
        except ValueError:
            relations.append(r)
            continue

        for bd in block_dims:
            # Avoid e.g. `x > yb`
            if any(i._depth > bd._depth for i in r[:n] if i.is_Block) or \
               any(bd._depth < i._depth for i in r[n+1:] if i.is_Block):
                continue

            relations.append(tuple(bd if i is d else i for i in r))

    # 3: Make sure BlockDimensions at same depth stick next to each other
    # E.g., `(t, xbb, ybb, xb, yb, x, y)`, and NOT e.g. `(t, xbb, xb, x, ybb, ...)`
    # NOTE: this is perfectly legal since:
    # TILABLE => (perfect nest & PARALLEL) => interchangeable
    for i in ispace.itdimensions:
        if not i.is_Block:
            continue
        for bd in block_dims:
            if bd._depth < i._depth:
                relations.append((bd, i))

    intervals = IntervalGroup(intervals, relations=relations)

    sub_iterators = dict(ispace.sub_iterators)
    sub_iterators.pop(d, None)
    sub_iterators.update({block_dims[-1]: ispace.sub_iterators.get(d, [])})
    sub_iterators.update({bd: () for bd in block_dims[:-1]})

    directions = dict(ispace.directions)
    directions.pop(d)
    directions.update({bd: ispace.directions[d] for bd in block_dims})

    return IterationSpace(intervals, sub_iterators, directions)


def skewing(clusters, options):
    """
    This pass helps to skew accesses and loop bounds as well as perform loop interchange
    towards wavefront temporal blocking
    Parameters
    ----------
    clusters : tuple of Clusters
        Input Clusters, subject of the optimization pass.
    options : dict
        The optimization options.
        * `skewinner` (boolean, False): enable/disable loop skewing along the
           innermost loop.
    """
    processed = clusters
    if options['blocktime']:
        processed = TBlocking(options).process(processed)

    processed = Skewing(options).process(processed)
    processed = RelaxSkewed(options).process(processed)

    return processed


class Skewing(Queue):

    """
    Construct a new sequence of clusters with skewed expressions and iteration spaces.

    Notes
    -----
    This transformation is applying loop skewing to derive the
    wavefront method of execution of nested loops. Loop skewing is
    a simple transformation of loop bounds and is combined with loop
    interchanging to generate the wavefront [1]_.

    .. [1] Wolfe, Michael. "Loops skewing: The wavefront method revisited."
    International Journal of Parallel Programming 15.4 (1986): 279-293.

    Examples:

    .. code-block:: python

        for i = 2, n-1
            for j = 2, m-1
                a[i,j] = (a[a-1,j] + a[i,j-1] + a[i+1,j] + a[i,j+1]) / 4

    to

    .. code-block:: python

        for i = 2, n-1
            for j = 2+i, m-1+i
                a[i,j-i] = (a[a-1,j-i] + a[i,j-1-i] + a[i+1,j-i] + a[i,j+1-i]) / 4

    """

    template = "%s%d_blk%s"

    def __init__(self, options):
        self.skewinner = bool(options['blockinner'])
        self.levels = options['blocklevels']

        super().__init__()

    def callback(self, clusters, prefix):
        if not prefix:
            return clusters

        d = prefix[-1].dim
        processed = []
        for c in clusters:
            if SKEWABLE not in c.properties[d]:
                return clusters

            if d is c.ispace[-1].dim and not self.skewinner:
                return clusters

            skew_dims = [i.dim for i in c.ispace if SEQUENTIAL in c.properties[i.dim]]
            if len(skew_dims) > 2:
                return clusters
            skew_dim = skew_dims[-1]

            # Since we are here, prefix is skewable and nested under a
            # SEQUENTIAL loop.

            skewlevel = 1
            intervals = []
            for i in c.ispace:
                if i.dim is d:
                    # If time is blocked skew at skewlevel + 1
                    cond1 = len(skew_dims) == 2 and d._depth == skewlevel + 1
                    # If time is blocked skew at level == 0 (e.g. subdims)
                    cond3 = len(skew_dims) == 2 and d._depth == 0
                    # If time is not blocked skew at level <=1
                    cond2 = len(skew_dims) == 1 and d._depth <= skewlevel

                    if cond1:
                        intervals.append(Interval(d, i.lower, i.upper))
                    elif cond2 or cond3:
                        intervals.append(Interval(d, skew_dim, skew_dim))
                    else:
                        intervals.append(i)
                else:
                    intervals.append(i)

            intervals = IntervalGroup(intervals, relations=c.ispace.relations)
            ispace = IterationSpace(intervals, c.ispace.sub_iterators,
                                    c.ispace.directions)

            exprs = xreplace_indices(c.exprs, {d: d - skew_dim})
            processed.append(c.rebuild(exprs=exprs, ispace=ispace,
                                       properties=c.properties))

        return processed


# Utils


def attach_property(clusters, d, p):
    """
    Attach `p` to `c.properties[d]` for each `c` in `clusters`.
    """
    processed = []
    for c in clusters:
        properties = dict(c.properties)
        properties[d] = set(properties[d]) | {p}
        processed.append(c.rebuild(properties=properties))

    return processed


class TBlocking(Queue):

    template = "%s%d_blk%s"

    def __init__(self, options):
        self.nblocked = Counter()
        super(TBlocking, self).__init__()

    def callback(self, clusters, prefix):
        if not prefix:
            return clusters

        d = prefix[-1].dim

        processed = []

        for c in clusters:
            sf = get_skewing_factor(c)
            if d.is_Time:
                name = self.template % (d.name, self.nblocked[d], '%d')
                block_dims = create_block_dims(name, d, 1, sf=sf)

                ispace = decompose(c.ispace, d, block_dims)
                # Use the innermost IncrDimension in place of `d`
                exprs = [uxreplace(e, {d: block_dims[-1]}) for e in c.exprs]

                # The new sub_iterators (skewing factor)
                sub_iterators = dict(ispace.sub_iterators)
                sub_iters = []
                for j in sub_iterators[block_dims[-1]]:
                    if sf > 1 and j.is_Modulo:
                        nom = INT(block_dims[-1]/sf)
                        denom = (block_dims[-1].root.symbolic_max -
                                 block_dims[-1].root.symbolic_min + 1)
                        sub_iters.append(j.func(offset=Mod(nom, denom) + j.offset - d))
                    else:
                        sub_iters.append(j)
                sub_iterators.update({block_dims[-1]: tuple(sub_iters)})

                # Should use corect intervals
                ispace = IterationSpace(ispace.intervals, sub_iterators,
                                        ispace.directions)

                # The new Cluster properties
                properties = dict(c.properties)
                properties.pop(d)
                properties.update({bd: c.properties[d] for bd in block_dims})

                processed.append(c.rebuild(exprs=exprs, ispace=ispace,
                                           properties=properties))

            elif d._depth == 1:  # Interchanged non-Time loops are not PARALLEL anymore
                properties = dict(c.properties)
                properties.update({d: c.properties[d] - {PARALLEL}})
                processed.append(c.rebuild(properties=properties))
            else:
                processed.append(c)
        return processed


class RelaxSkewed(Queue):

    def __init__(self, options):
        self.nblocked = Counter()
        super(RelaxSkewed, self).__init__()

    def callback(self, clusters, prefix):
        if not prefix:
            return clusters

        d = prefix[-1].dim

        # Rule out time dim and not is_Incr
        if d.is_Time or not d.is_Incr:
            return clusters

        processed = []
        for c in clusters:
            skew_dims = [i.dim for i in c.ispace if SEQUENTIAL in c.properties[i.dim]]

            if len(skew_dims) == 1:
                processed.append(c)
                continue

            family_dims = [j.dim for j in c.ispace if j.dim.root is d.root]

            # IMPORTANT: we only process the head of the family in this Queue
            if d is not family_dims[0]:
                processed.append(c)
                continue

            sf = get_skewing_factor(c)
            skew_dim = skew_dims[-1]
            intervals = []
            mapper = {}
            for i in c.ispace:
                if i.dim in family_dims:
                    if i.dim._depth == 1:
                        offset = sf*(skew_dim.root.symbolic_max -
                                     skew_dim.root.symbolic_min)
                        rmax = i.dim.symbolic_max + offset
                        sd = i.dim.func(rmax=rmax)
                        intervals.append(Interval(sd, i.lower, i.upper))
                        mapper.update({i.dim: sd})
                    elif i.dim._depth == 2:
                        rmin = evalrel(max, [i.dim.symbolic_min,
                                       i.dim.root.symbolic_min + skew_dim])
                        rmax = i.dim.symbolic_rmax.xreplace({i.dim.root.symbolic_max:
                                                            i.dim.root.symbolic_max +
                                                            skew_dim})
                        sd2 = i.dim.func(parent=sd, rmin=rmin, rmax=rmax)
                        intervals.append(Interval(sd2, i.lower, i.upper))
                        mapper.update({i.dim: sd2})
                    elif i.dim._depth > 2:
                        res = evalrel(min, [sd2.symbolic_rmax, i.dim.symbolic_max])
                        rmax = i.dim.symbolic_rmax.xreplace({i.dim.symbolic_rmax: res})
                        sd3 = i.dim.func(parent=sd2, rmax=rmax)
                        intervals.append(Interval(sd3, i.lower, i.upper))
                        mapper.update({i.dim: sd3})
                        sd2 = sd3
                    else:
                        intervals.append(i)
                else:
                    intervals.append(i)

            # After rebase relations are now empty , to fix
            # Update `relations` with the newly created `Dimension`s
            relations = []
            for r in c.ispace.relations:
                if any(f in r for f in family_dims) and mapper:
                    rl = as_list(r)
                    newr = [j.xreplace(mapper) for j in rl]
                    relations.append(as_tuple(newr))
                else:
                    relations.append(r)

            # Sanity check
            assert len(relations) == len(c.ispace.relations)
            # Build new intervals
            intervals = IntervalGroup(intervals, relations=relations)

            # Update `sub_iterators`, `directions`, `properties`, `expressions`
            sub_iterators = dict(c.ispace.sub_iterators)
            directions = dict(c.ispace.directions)
            properties = dict(c.properties)

            for f in family_dims:
                sub_iterators.pop(f, None)
                sub_iterators.update({mapper[f]: c.ispace.sub_iterators.get(f, [])})
                directions.pop(f)
                directions.update({mapper[f]: c.ispace.directions[f]})
                properties.pop(f)
                properties.update({mapper[f]: c.properties[f]})
                exprs = xreplace_indices(c.exprs, {f: mapper[f]})

            # Build the new `IterationSpace`
            ispace = IterationSpace(intervals, sub_iterators, directions)
            processed.append(c.rebuild(exprs=exprs, ispace=ispace,
                                       properties=properties))

        return processed


def get_skewing_factor(cluster):
    '''
    Returns the skewing factor needed to skew a cluster of expressions. Skewing factor is
    equal to half the maximum of the functions' space orders and helps to preserve valid
    data dependencies while skewing.

    Parameters
    ----------
    cluster: Cluster
        Input Cluster, subject of the computation
    '''
    functions = {i.function for i in retrieve_indexed(cluster.exprs)}
    try:
        sf = int(max([i.space_order for i in functions])/2)
    except AttributeError:
        sf = 1
    return (sf if sf else 1)
