# Code repository for thesis Chapter 3

# Prerequisites:
  - A working devito installation https://www.devitoproject.org/devito/download.html
  
# To install in Python virual environment:
  - Download and extract the code from https://zenodo.org/record/6962402
  - Create a new python virtual environment
  - `cd georgebisbas-devito-...`
  - `pip3 install -e .[extras]`
  - Set env variables and pinning.

```bash
export DEVITO_LANGUAGE=openmp
export OMP_PROC_BIND=close
export DEVITO_LOGGING=DEBUG
export DEVITO_ARCH=gcc
```

Run :
```bash
DEVITO_JIT_BACKDOOR=0 python3 examples/seismic/acoustic/demo_temporal_sources.py -so 4
```

This will generate code for a space order 4 acoustic devito kernel with standard space blocking 
and another temporal blocking kernel that we will modify.
You will notice difference between `norm(usol)` and `norm(uref)`.
This is expected as the generated kernel needs to be manually modified.
The generated log will end executing a kernel under `===Temporal blocking================================`

```
===Temporal blocking======================================
Allocating memory for n(1,)
Allocating memory for usol(3, 236, 236, 236)
gcc -O3 -g -fPIC -Wall -std=c99 -march=native -Wno-unused-result -Wno-unused-variable -Wno-unused-but-set-variable -ffast-math -shared -fopenmp /tmp/devito-jitcache-uid1000/xxx-hash-xxx.c -lm -o /tmp/devito-jitcache-uid1000/xxx-hash-xxx.so
Operator `Kernel` jit-compiled `/tmp/devito-jitcache-uid1000/xxx-hash-xxx.c` in 0.48 s with `GNUCompiler`
Operator `Kernel` run in 0.49 s
```

Within the kernels folders under the following folders we have stored the manually edited kernels with temporal blocking:

```bash
devito/examples/seismic/acoustic/kernels/
devito/examples/seismic/elastic/kernels/
devito/examples/seismic/tti/kernels/
```

For example:
Copy the space order 4 kernel from `devito/examples/seismic/acoustic/kernels/hash` to the `hash-xxx.c`
`cp kernels/8393*...xxx*.c /tmp/devito-jitcache-uid1000/8393*.c`

Then try `DEVITO_JIT_BACKDOOR=1 python3 examples/seismic/acoustic/demo_temporal_sources.py -so 4`
to re-run. Now the norms should match and we are ready to run the experiemnts.
If not contact me ASAP :-).

Use arguments `-d nx ny nx` , `-tn timesteps` to pass as arguments domain size and number of timesteps.
e.g.:
`DEVITO_JIT_BACKDOOR=1 python3 examples/seismic/acoustic/demo_temporal_sources.py -d 512 512 512 --tn 512 -so 8`

The available kernels are for space orders 4, 8, 12.


## Experiments for Chapter3
```bash
export DEVITO_ARCH=gcc
export DEVITO_LOGGING=DEBUG
export DEVITO_LANGUAGE=openmp
```
### Acoustic kernel
```bash
DEVITO_JIT_BACKDOOR=0 python3.8 examples/seismic/acoustic/acoustic_tune_model.py -so 4 -d 512 512 512 --tn 512
will fail... follow instructions above to copy the correct kernel

# After the correct kernel has been copied:
DEVITO_JIT_BACKDOOR=1 python3.8 examples/seismic/acoustic/acoustic_tune_model.py -so 4 -d 512 512 512 --tn 512
DEVITO_JIT_BACKDOOR=1 python3.8 examples/seismic/acoustic/acoustic_tune_model.py -so 4 -d 512 512 512 --tn 512
DEVITO_JIT_BACKDOOR=1 python3.8 examples/seismic/acoustic/acoustic_tune_model.py -so 4 -d 512 512 512 --tn 512
```

### Elastic kernel

```bash
DEVITO_JIT_BACKDOOR=0 python3.8 examples/seismic/elastic/elastic_tune_model.py -so 4 -d 512 512 512 --tn 512
will fail... follow instructions above to copy the correct kernel

# After the correct kernel has been copied:
DEVITO_JIT_BACKDOOR=1 python3.8 examples/seismic/elastic/elastic_tune_model.py -so 4 -d 512 512 512 --tn 512
DEVITO_JIT_BACKDOOR=1 python3.8 examples/seismic/elastic/elastic_tune_model.py -so 4 -d 512 512 512 --tn 512
DEVITO_JIT_BACKDOOR=1 python3.8 examples/seismic/elastic/elastic_tune_model.py -so 4 -d 512 512 512 --tn 512
```

### TTI kernel

```bash
DEVITO_JIT_BACKDOOR=0 python3.8 examples/seismic/tti/tti_tune_model.py -so 4 -d 512 512 512 --tn 512
will fail... follow instructions above to copy the correct kernel

# After the correct kernel has been copied:
DEVITO_JIT_BACKDOOR=1 python3.8 examples/seismic/tti/tti_tune_model.py -so 4 -d 512 512 512 --tn 512
DEVITO_JIT_BACKDOOR=1 python3.8 examples/seismic/tti/tti_tune_model.py -so 4 -d 512 512 512 --tn 512
DEVITO_JIT_BACKDOOR=1 python3.8 examples/seismic/tti/tti_tune_model.py -so 4 -d 512 512 512 --tn 512
```


# Tuning Devito
Run Devito with `DEVITO_LOGGING=aggressive` so as to ensure that one of the best space-blocking configurations are selected.

# Tuning time-tiled kernel
In order to manually tune the Devito time-tiled kernel one should use an editor and jit-backdoor.
You can and should manually play around tile and block size values.
```
  int xb_size = 32;
  int yb_size = 32; // to fix as 8/16 etc
  int num_threads = 8;
  int x0_blk0_size = 8;
  int y0_blk0_size = 8;
```
You should also change accordingly the number of threads manually.
According to experiments x0_blk0_size, y0_blk0_size are best at 8
while xb_size, yb_size are nice in {32, 48, 64} depending on the platform.

`xb_size, yb_size == tiles`
`x0_blk0_size, y0_blk0_size == blocks`

When on Skylake: you may also need to change SIMD parallelism from 32 to 64.

As of YASK:
`From YASK: Although the terms "block" and "tile" are often used interchangeably, in
this section, we [arbitrarily] use the term "block" for spatial-only grouping
and "tile" when multiple temporal updates are allowed.`

Let me know your findings and your performance results. https://opesci-slackin.now.sh/ at #time-tiling

Depending on platform, I would expect speed-ups along the following lines:
![Perf results](https://github.com/devitocodes/devito/blob/timetiling_on_cd/examples/seismic/acoustic/kernels/temporal_performance.pdf)
