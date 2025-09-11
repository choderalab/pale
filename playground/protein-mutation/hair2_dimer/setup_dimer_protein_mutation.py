#!/usr/bin/env python
# coding: utf-8

import os
import re
from gufe import ChemicalSystem, ProteinComponent, SolventComponent, LigandAtomMapping
from feflow.protocols import NonEquilibriumCyclingProtocol

# TODO: Translate this into a shell(?) script to be run with `openfe` cli.

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

## UTILITY FUNCTIONS AND OBJECTS

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


# Protein components
structures_rootdir = "."
initial_structure = ProteinComponent.from_pdb_file(
    f"{structures_rootdir}/her2-for-mapping-p61r.pdb"
)
final_structure = ProteinComponent.from_pdb_file(f"{structures_rootdir}/mutated_dimer_P61R.pdb")
# Solvent component
solvent_comp = SolventComponent()
# Create chemical systems
initial_state = ChemicalSystem(components={"protein": initial_structure, "solvent": solvent_comp})
end_state = ChemicalSystem(components={"protein": final_structure, "solvent": solvent_comp})
# Mappings
mappings_rootdir = "./mappings"
input_file = f"{mappings_rootdir}/dimer_P61R.json"
with open(input_file) as in_file:
    mapping = LigandAtomMapping.from_json(in_file)

# Settings
settings = NonEquilibriumCyclingProtocol.default_settings()  # getting default settings
# changing default settings as needed
settings.integrator_settings.equilibrium_steps = 250
settings.integrator_settings.nonequilibrium_steps = 250
settings.alchemical_settings.explicit_charge_correction = True
# Better charge method with Openeye
settings.partial_charge_settings.partial_charge_method = "am1bccelf10"
settings.partial_charge_settings.off_toolkit_backend = "openeye"

# Define protocol/simulation
protocols_rootdir = "./protocol_dags_apo"
os.makedirs(protocols_rootdir, exist_ok=True)
protocol = NonEquilibriumCyclingProtocol(settings=settings)
protocol_dag = protocol.create(
    stateA=initial_state, stateB=end_state, mapping=mapping, name=f"dimer_P61R"
)
protocol_dag_path = f"{protocols_rootdir}/NEq_cycling_dimer_P61R.json"
protocol_dag.to_json(protocol_dag_path)
