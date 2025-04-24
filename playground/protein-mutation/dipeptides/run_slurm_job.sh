#!/bin/bash

#SBATCH --partition gpu
#SBATCH --ntasks=1
#SBATCH --array=1-20
#SBATCH --time=12:00:00
#SBATCH --job-name=neq-dipeptides-capped
#SBATCH --output=stdout/%A_%a_%x_%N.out
#SBATCH --error=stderr/%A_%a_%x_%N.err
#SBATCH --mem-per-cpu=8G
#SBATCH --gpus-per-task=1

# Print hostname
echo $SLURM_SUBMIT_HOST

# Source bashrc (useful for working conda/mamba)
source ${HOME}/.bashrc
OPENMM_CPU_THREADS=1

# Activate environment
conda activate feflow-openfe-env

# Report node in use
hostname

# Open eye license activation/env
export OE_LICENSE=${HOME}/.OpenEye/oe_license.txt

# Report CUDA info
env | sort | grep 'CUDA'

# set loglevel debug
LOGLEVEL=DEBUG

# launching a benchmark pair (target, edge) per job (0-based thus substract 1)
python run_protocol_dag.py --protocol-dags-dir ./protocol_dags_capped --index $(( ${SLURM_ARRAY_TASK_ID} - 1 ))

