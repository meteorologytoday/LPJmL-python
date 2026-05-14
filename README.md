# LPJmL-python

A standalone Python/C extension that drives LPJmL through its FMS coupler interface,
mimicking the POEM atmosphere–land coupling loop without needing the full POEM model.
Useful for testing and debugging LPJmL in isolation.

## Overview

The extension (`LPJmLCoupler`) wraps LPJmL's C driver functions so they can be called
from Python with NumPy arrays.  A typical simulation loop looks like:

```
lpjml_init → [ lpjml_update × (years × months × days × hours) ] → lpjml_end
```

All per-cell fields (forcings and outputs) are passed as 1-D NumPy arrays of length `ncells`.

## Repository layout

```
src/
  python_interface.c   # Python/C extension wrapping the LPJmL driver
  lpj_standalone.c     # Reference copy of the LPJmL driver (do not edit)
testrun/
  main.py              # Example driver script
  lpjml_fms.cjson      # Full-run config  (1840–1841, reads from restart)
  lpjml_fms_spinup.cjson   # Spinup config
  lpjml_fms_test100.cjson  # Short test config
  run.slurm.sh         # SLURM submission script
setup.py               # Build configuration
compile.sh             # Convenience build script
```

## Prerequisites

- Intel MPI compiler (`mpiicc`)
- Pre-built LPJmL installation at the path set in `setup.py` (`LPJROOT`)
- Python environment with NumPy (available via the POEM module system)

## Building the extension

Load the environment first, then build in-place:

```bash
loadpoem artificial_wet_soil_moisture
bash compile.sh
```

This runs `python setup.py build_ext --inplace` and produces
`LPJmLCoupler.cpython-*.so` in the project root.

## Running `testrun/main.py`

### Directly (login node or interactive session)

```bash
cd testrun
export PYTHONPATH=/p/projects/poem/tienyiao/projects/LPJmL-python:$PYTHONPATH
python main.py lpjml_fms.cjson
```

The config file argument is optional; it defaults to `lpjml_fms.cjson`.

### Via SLURM

```bash
cd testrun
sbatch run.slurm.sh [config.cjson]
```

`run.slurm.sh` defaults to `lpjml_fms_spinup.cjson` when no config is given.
Logs are written to `sbatch.64.<jobid>.{out,err}`.

### With MPI (multiple ranks)

```bash
mpirun -n 4 python main.py lpjml_fms.cjson
```

## What `main.py` does

1. **Init** — calls `lpjml_init(config, dt_fast=3600)`, which returns grid metadata
   (`ncells`, `firstyear`, `lastyear`, `nspinup`).
2. **Allocate arrays** — creates NumPy arrays for atmospheric forcings
   (precipitation, temperature, shortwave/longwave radiation, wind, CO₂, …)
   and output fields (carbon flux, runoff, albedo, surface temperature, …).
3. **Time loop** — iterates over years → months → days → hours, calling
   `lpjml_update(...)` once per hour with the current forcings.
   The previous timestep's canopy air humidity (`q_ca`) is fed back as input.
4. **Finalise** — calls `lpjml_end()`.

## Available configs

| File | Spinup years | Sim years | Notes |
|------|-------------|-----------|-------|
| `lpjml_fms_spinup.cjson` | set in file | — | Spinup run |
| `lpjml_fms_test100.cjson` | — | short | Quick test |
| `lpjml_fms.cjson` | 0 | 1840–1841 | Reads from `restart/restart_POEM_2022.lpj` |
