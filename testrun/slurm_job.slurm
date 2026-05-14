#! /bin/sh
##SBATCH --account=next
#SBATCH --job-name=LPJ-Forcing-interp
#
#SBATCH --tasks 64
#SBATCH --nodes 1
#
# Acceptable time formats include "minutes", "minutes:seconds", "hours:minutes:seconds",
#                   "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"
#SBATCH --qos=short
#SBATCH --time=23:59:00 
#
#SBATCH --account=poem
#SBATCH --output=sbatch.64.%j.out
#SBATCH --error=sbatch.64.%j.err
#SBATCH --mail-type=TIME_LIMIT


source /home/tienyiao/.bashrc_POEM 3
export PYTHONPATH=/p/projects/poem/tienyiao/projects/LPJmL-python:$PYTHONPATH

PYTHON=/home/tienyiao/miniconda3/envs/rainbow/bin/python

NJOBS=${SLURM_NTASKS:-20}

mpirun -n $NJOBS $PYTHON $SLURM_SUBMIT_DIR/main.py $SLURM_SUBMIT_DIR/lpjml_fms.cjson
