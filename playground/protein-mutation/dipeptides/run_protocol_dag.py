#!/usr/bin/env python
# coding: utf-8

import argparse
import pathlib
import re
from feflow.protocols import NonEquilibriumCyclingProtocol
from gufe.protocols import execute_DAG
# Script to execute openfe/feflow protocol dags given their serialized json file

# Basic CLI
parser = argparse.ArgumentParser(description="Run protocol dag given index")
parser.add_argument(
    "--protocol-dags-dir",
    type=str,
    help="base directory where serialized protocol dags are located",
)
parser.add_argument("--index", type=int, help="index of protocol dag to file")

args = parser.parse_args()

# Helper objects/functions
# Mutation specification
mutations_to_perform = [
    "Y2F",
    "Y2A",
    "W2F",
    "T2A",
    "E2A",
    "D2A",
    "K2A",
    "R2A",
    "R2Q",
    "H2A",
    "F2Y",
    "A2Y",
    "F2W",
    "A2T",
    "A2E",
    "A2D",
    "A2K",
    "A2R",
    "Q2R",
    "A2H",
]

# Create dictionary for AA code translation
aa_three_to_one_code = {
    "ALA": "A",
    "GLY": "G",
    "ILE": "I",
    "LEU": "L",
    "PRO": "P",
    "VAL": "V",
    "PHE": "F",
    "TRP": "W",
    "TYR": "Y",
    "ASP": "D",
    "GLU": "E",
    "ARG": "R",
    "HIS": "H",
    "LYS": "K",
    "SER": "S",
    "THR": "T",
    "CYS": "C",
    "MET": "M",
    "ASN": "N",
    "GLN": "Q",
}
aa_one_to_three_code = {value: key for key, value in aa_three_to_one_code.items()}


def parse_mutation_spec(mutation_spec: str):
    """
    Parse a mutation specification string into its components: initial residue,
    residue number, and final residue.

    The input is expected to follow the format:
    "<initial_residue><residue_number><final_residue>", where both residues
    can be one- or three-letter uppercase or lowercase codes. The function
    automatically capitalizes the input before parsing.

    Examples of valid formats:
    - "Y2A"
    - "tyr23ala"
    - "Ace123NME"

    Parameters
    ----------
    mutation_spec : str
        A string representing a mutation (e.g., "Y2A", "TYR23ALA", "ace123nme").

    Returns
    -------
    tuple of str
        A tuple (initial_residue, residue_number, final_residue), where:
            - initial_residue is a 1–3 character uppercase string,
            - residue_number is a string of digits,
            - final_residue is a 1–3 character uppercase string.

    Raises
    ------
    ValueError
        If the mutation_spec does not match the expected pattern.
    """
    mutation_string = mutation_spec.upper()
    pattern = r"([A-Z]{1,3})(\d+)([A-Z]{1,3})"
    match = re.search(pattern, mutation_string)
    if match:
        initial_aa, res_number, final_aa = match.groups()
        return initial_aa, res_number, final_aa
    raise ValueError(f"Invalid mutation specification: {mutation_spec}")


# Main
protocol_dags_dir = args.protocol_dags_dir
protocol_dag_index = args.index
# Enumerating mutations -- this is just useful to be used by slurm jobs array ids
mutations_enum = dict(enumerate(mutations_to_perform))
print(mutations_enum)
initial_aa, _, final_aa = parse_mutation_spec(mutations_enum[protocol_dag_index])
# converting to three letter code -- TODO: Make this a util callable
initial_aa = aa_one_to_three_code[initial_aa]
final_aa = aa_one_to_three_code[final_aa]
# Execute protocol dag
# read from serialized protocol dag
protocol_dag_path = f"./{protocol_dags_dir}/NEq_cycling_{initial_aa}_to_{final_aa}.json"
protocol_dag_deserialized = NonEquilibriumCyclingProtocol.from_json(protocol_dag_path)
# Execute dag
# TODO: Make the results_path part of the CLI args
results_path = pathlib.Path("./results_capped")
results_path.mkdir(exist_ok=True)
print(f"Executing protocol dag for {initial_aa} to {final_aa} mutation.")
protocol_result_dag = execute_DAG(
    protocol_dag_deserialized,
    keep_shared=True,
    shared_basedir=results_path,
    scratch_basedir=results_path,
)
protocol_result_dag.to_json(f"./{results_path}/results_NEq_cycling_{initial_aa}_to_{final_aa}.json")
