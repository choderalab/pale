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

args = parser.parse_args()


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
# Execute protocol dag
# read from serialized protocol dag
protocol_dag_path = f"./{protocol_dags_dir}/NEq_cycling_dimer_P61R.json"
protocol_dag_deserialized = NonEquilibriumCyclingProtocol.from_json(protocol_dag_path)
# Execute dag
# TODO: Make the results_path part of the CLI args
results_path = pathlib.Path("./results_capped")
results_path.mkdir(exist_ok=True)
print(f"Executing protocol dag.")
protocol_result_dag = execute_DAG(
    protocol_dag_deserialized,
    keep_shared=True,
    shared_basedir=results_path,
    scratch_basedir=results_path,
)
protocol_result_dag.to_json(f"./{results_path}/results_NEq_cycling_dimer_P61R.json")
