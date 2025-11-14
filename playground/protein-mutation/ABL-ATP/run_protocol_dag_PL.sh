#!/bin/bash
#SBATCH --partition=componc_gpu_batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=16GB
#SBATCH --time=2:00:00
#SBATCH --gpus=1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=pathilm@mskcc.org

SCRIPT_DIRECTORY_PATH=/data1/choderaj/pathilm/FEC/kinase_ATP
ENV_NAME=pale-espaloma

source ~/.bashrc
conda activate $ENV_NAME

if [ -n "${INPUT_SDF+x}" ]; then
    echo "HOLO: Running script with ligand included."
    python run_protocol_dag_PL.py -p $INPUT_PDB -l $INPUT_SDF -m $MUTATION --leg $LEG --num_cycles $NUM_CYCLES
else
    echo "APO: Running script without ligand."
    python run_protocol_dag_PL.py -p $INPUT_PDB -m $MUTATION --leg $LEG --num_cycles $NUM_CYCLES
fi