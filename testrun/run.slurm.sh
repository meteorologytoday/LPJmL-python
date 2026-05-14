#!/bin/bash
#SBATCH --job-name=LPJ-Forcing-interp
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --qos=short
#SBATCH --time=23:59:00
#SBATCH --account=poem
#SBATCH --output=sbatch.64.%j.out
#SBATCH --error=sbatch.64.%j.err
#SBATCH --mail-type=TIME_LIMIT

# Usage:
#   As a SLURM job:  sbatch run.sh [config.cjson]
#   On a login node: bash run.sh [config.cjson]

source /home/tienyiao/.bashrc_POEM 3
export PYTHONPATH=/p/projects/poem/tienyiao/projects/LPJmL-python:$PYTHONPATH

PYTHON=/p/system/packages/tools/miniforge/25.3.0-3/bin/python3

# Config file: first argument, or default
CONFIG=${1:-lpjml_fms_spinup.cjson}

# When submitted via sbatch, resolve paths relative to submission directory.
# When run locally, use the directory of this script.
if [ -n "$SLURM_SUBMIT_DIR" ]; then
    RUNDIR=$SLURM_SUBMIT_DIR
else
    RUNDIR=$(cd "$(dirname "$0")" && pwd)
fi

NJOBS=${SLURM_NTASKS:-1}

mpirun -n $NJOBS $PYTHON $RUNDIR/main.py $RUNDIR/$CONFIG
