#!/bin/bash
#SBATCH --job-name="abl_rbfe"
#SBATCH --partition=cpushort
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=2GB
#SBATCH --time=5:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=pathilm@mskcc.org

SCRIPT_DIRECTORY_PATH=/data1/choderaj/pathilm/FEC/kinase_ATP
ENV_NAME=pale-espaloma

source ~/.bashrc
conda activate $ENV_NAME

PROTEIN=/data1/choderaj/pathilm/FEC/kinase_ATP/abl/structures/kinase_final.pdb
LIGAND=/data1/choderaj/pathilm/FEC/kinase_ATP/abl/structures/ATP_final.sdf

MUTATIONS="T315I G250E Y253F F317L N368S"
NUM_CYCLES=20


for MUTATION in $MUTATIONS; do
    sbatch \
    --time=20:00:00 \
    --export=INPUT_PDB=$PROTEIN,MUTATION=$MUTATION,LEG="apo",NUM_CYCLES=$NUM_CYCLES \
    --output=/data1/choderaj/pathilm/FEC/kinase_ATP/abl/$MUTATION-apo-%a.log \
    --error=/data1/choderaj/pathilm/FEC/kinase_ATP/abl/$MUTATION-apo-%a.stderr \
    run_protocol_dag_PL.sh

    sbatch \
    --time=20:00:00 \
    --export=INPUT_PDB=$PROTEIN,INPUT_SDF=$LIGAND,MUTATION=$MUTATION,LEG="holo",NUM_CYCLES=$NUM_CYCLES \
    --output=/data1/choderaj/pathilm/FEC/kinase_ATP/abl/$MUTATION-holo-%a.log \
    --error=/data1/choderaj/pathilm/FEC/kinase_ATP/abl/$MUTATION-holo-%a.stderr \
    run_protocol_dag_PL.sh
done
