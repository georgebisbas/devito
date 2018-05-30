from __future__ import absolute_import

from collections import OrderedDict
from itertools import combinations
from functools import reduce
from operator import mul
import resource

from devito.ir.iet import Iteration, FindNodes, FindSymbols
from devito.logger import info, info_at
from devito.parameters import configuration

__all__ = ['autotune']


def autotune(operator, arguments, tunable):
    """
    Acting as a high-order function, take as input an operator and a list of
    operator arguments to perform empirical autotuning. Some of the operator
    arguments are marked as tunable.
    """
    at_arguments = arguments.copy()

    # User-provided output data must not be altered
    output = [i.name for i in operator.output]
    for k, v in arguments.items():
        if k in output:
            at_arguments[k] = v.copy()

    iterations = FindNodes(Iteration).visit(operator.body)
    dim_mapper = {i.dim.name: i.dim for i in iterations}

    # Shrink the iteration space of time-stepping dimension so that auto-tuner
    # runs will finish quickly
    steppers = [i for i in iterations if i.dim.is_Time or i.dim.is_Stepping]
    if len(steppers) == 0:
        timesteps = 1
    elif len(steppers) == 1:
        stepper = steppers[0]
        start = stepper.dim.rtargs.start.default_value
        timesteps = stepper.extent(start=start, finish=options['at_squeezer'])
        if timesteps < 0:
            timesteps = options['at_squeezer'] - timesteps + 1
            info_at("Adjusted auto-tuning timestep to %d" % timesteps)
        at_arguments[stepper.dim.symbolic_start.name] = start
        at_arguments[stepper.dim.symbolic_end.name] = timesteps
        if stepper.dim.is_Stepping:
            at_arguments[stepper.dim.parent.symbolic_start.name] = start
            at_arguments[stepper.dim.parent.symbolic_end.name] = timesteps
    else:
        info_at("Couldn't understand loop structure, giving up auto-tuning")
        return arguments

    # Attempted block sizes ...
    mapper = OrderedDict([(i.argument.symbolic_size.name, i) for i in tunable])
    time_dim = None
    for i, d in mapper.items():
        if d.original_dim.is_Time:
            time_dim = i

    # ... Defaults (basic mode)
    blocksizes = [OrderedDict([(i, v) for i in mapper if not mapper[i].original_dim.is_Time]) for v in options['at_blocksize']]  # cubes
    # ... Always try the entire iteration space (degenerate block)
    datashape = [at_arguments[mapper[i].original_dim.symbolic_end.name] -
                 at_arguments[mapper[i].original_dim.symbolic_start.name] for i in mapper]
    blocksizes.append(OrderedDict([(i, mapper[i].iteration.extent(0, j))
                      for i, j in zip(mapper, datashape)]))  # degenerate block
    # ... More attempts if auto-tuning in aggressive mode
    if configuration.core['autotuning'] == 'aggressive':
        last_dim = None
        innermost = iterations[-1].dim
        for k, v in mapper.items():
            if v.original_dim == innermost:
                last_dim = (k, blocksizes[-1][k])

        blocksizes = more_heuristic_attempts(blocksizes)

        if last_dim:
            info_at("Extending the innermost dimension, %s <%s>" % (last_dim[0], last_dim[1]))
            intermediate_blocks = [OrderedDict([(i, v) for i in mapper if not (mapper[i].original_dim.is_Time or mapper[i].original_dim == innermost)])
                                   for v in options['at_blocksize']]
            intermediate_blocks = more_heuristic_attempts(intermediate_blocks)
            blocksizes += cross_time_tiles(intermediate_blocks, last_dim[0], [last_dim[1]])
            # TODO: don't extend this: run generator for 2 dims, then extend that

    if time_dim:
        blocksizes = cross_time_tiles(blocksizes, time_dim, [1, 2, 4, 8, 16])


    # How many temporaries are allocated on the stack?
    # Will drop block sizes that might lead to a stack overflow
    functions = FindSymbols('symbolics').visit(operator.body +
                                               operator.elemental_functions)
    stack_shapes = [i.shape for i in functions if i.is_Array and i._mem_stack]
    stack_space = sum(reduce(mul, i, 1) for i in stack_shapes)*operator.dtype().itemsize

    # Note: there is only a single loop over 'blocksize' because only
    # square blocks are tested
    timings = OrderedDict()
    fastest, timing = None, float("inf")
    unique = []

    for bs in blocksizes:
        if bs in unique:
            continue
        unique.append(bs)

        illegal = False
        for k, v in at_arguments.items():
            if k in bs:
                val = bs[k]
                start = at_arguments[mapper[k].original_dim.symbolic_start.name]
                end = at_arguments[mapper[k].original_dim.symbolic_end.name]
                if val <= mapper[k].iteration.extent(start, end):
                    at_arguments[k] = val
                else:
                    # Block size cannot be larger than actual dimension
                    illegal = True
                    break
        if illegal:
            continue

        # Make sure we remain within stack bounds, otherwise skip block size
        dim_sizes = {}
        for k, v in at_arguments.items():
            if k in bs:
                dim_sizes[mapper[k].argument.symbolic_size] = bs[k]
            elif k in dim_mapper:
                dim_sizes[dim_mapper[k].symbolic_size] = v
        try:
            bs_stack_space = stack_space.xreplace(dim_sizes)
        except AttributeError:
            bs_stack_space = stack_space
        try:
            if int(bs_stack_space) > options['at_stack_limit']:
                continue
        except TypeError:
            # We should never get here
            info_at("Couldn't determine stack size, skipping block size %s" % str(bs))
            continue

        # Use AT-specific profiler structs
        at_arguments[operator.profiler.varname] = operator.profiler.setup()

        operator.cfunction(*list(at_arguments.values()))
        elapsed = sum(operator.profiler.timings.values())
        timings[tuple(bs.items())] = elapsed
        if elapsed < timing:
            fastest = tuple(bs.items())
            timing = elapsed
        info_at("Block shape <%s> took %f (s) in %d time steps" %
                (','.join('%d' % i for i in bs.values()), elapsed, timesteps))

    try:
        # best = dict(min(timings, key=timings.get))
        best = dict(fastest)
        info("Auto-tuned block shape: %s; time: %f (s)" % (best, timing))
    except ValueError:
        info("Auto-tuning request, but couldn't find legal block sizes")
        return arguments

    # Build the new argument list
    tuned = OrderedDict()
    for k, v in arguments.items():
        tuned[k] = best[k] if k in mapper else v

    # Reset the profiling struct
    assert operator.profiler.varname in tuned
    tuned[operator.profiler.varname] = operator.profiler.setup()

    return tuned


def more_heuristic_attempts(blocksizes):
    # Ramp up to higher block sizes
    handle = OrderedDict([(i, options['at_blocksize'][-1]) for i in blocksizes[0]])
    # insert more cubes
    for i in range(3):
        new_bs = OrderedDict([(k, v*2) for k, v in handle.items()])
        blocksizes.insert(blocksizes.index(handle) + 1, new_bs)
        handle = new_bs

    handle = []
    # Extended shuffling for the smaller block sizes
    for bs in blocksizes[:4]:
        for i in blocksizes:
            handle.append(OrderedDict(list(bs.items())[:-1] + [list(i.items())[-1]]))
    # Some more shuffling for all block sizes
    for bs in list(blocksizes):
        ncombs = len(bs)  # dimensions to tile over
        for i in range(ncombs):
            for j in combinations(bs, i+1):
                item = [(k, bs[k]*2 if k in j else v) for k, v in bs.items()]
                handle.append(OrderedDict(item))

    return blocksizes + handle


def extend_dimension(blocksizes, dim, size):
    return blocksizes + [OrderedDict([(dim, size) if dim == d else (d, s) for d, s in bs.items()]) for bs in blocksizes]


def cross_time_tiles(blocksizes, dim, tiles):
    extended = []
    for bs in blocksizes:
        for tile in tiles:
            extended.append(OrderedDict([(dim, tile)] + list(bs.items())))

    return extended


options = {
    'at_squeezer': 17,
    'at_blocksize': sorted({8, 16, 24, 32, 40, 64, 128}),
    'at_stack_limit': resource.getrlimit(resource.RLIMIT_STACK)[0] / 4
}
"""Autotuning options."""
