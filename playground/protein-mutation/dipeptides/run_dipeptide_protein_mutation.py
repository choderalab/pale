#!/usr/bin/env python
# coding: utf-8

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


for mutation_spec in mutations_to_perform:
    initial_aa, resnumber, final_aa = parse_mutation_spec(mutation_spec=mutation_spec)

    # converting to three letter code
    initial_aa = aa_one_to_three_code[initial_aa]
    final_aa = aa_one_to_three_code[final_aa]

    # Protein components
    structures_rootdir = "./structures"
    initial_structure = ProteinComponent.from_pdb_file(
        f"{structures_rootdir}/{initial_aa}_capped.pdb"
    )
    final_structure = ProteinComponent.from_pdb_file(f"{structures_rootdir}/{final_aa}_capped.pdb")
    # Solvent component
    solvent_comp = SolventComponent()
    # Create chemical systems
    initial_state = ChemicalSystem(
        components={"protein": initial_structure, "solvent": solvent_comp}
    )
    end_state = ChemicalSystem(components={"protein": final_structure, "solvent": solvent_comp})
    # Mappings
    mappings_rootdir = "./mappings"
    input_file = f"{mappings_rootdir}/{initial_aa}_to_{final_aa}.json"
    with open(input_file) as in_file:
        mapping = LigandAtomMapping.from_json(in_file)

    # Settings
    settings = NonEquilibriumCyclingProtocol.default_settings()  # getting default settings
    # changing default settings as needed
    settings.integrator_settings.equilibrium_steps = 375000
    settings.integrator_settings.nonequilibrium_steps = 375000
    settings.alchemical_settings.explicit_charge_correction = True
    # Better charge method with Openeye
    settings.partial_charge_settings.partial_charge_method = "am1bccelf10"
    settings.partial_charge_settings.off_toolkit_backend = "openeye"

    # Define protocol/simulation
    protocols_rootdir = "./protocol_dags_capped"
    protocol = NonEquilibriumCyclingProtocol(settings=settings)
    protocol_dag = protocol.create(
        stateA=initial_state, stateB=end_state, mapping=mapping, name=f"{initial_aa}_to_{final_aa}"
    )
    protocol_dag_path = f"{protocols_rootdir}/NEq_cycling_{initial_aa}_to_{final_aa}.json"
    protocol_dag.to_json(protocol_dag_path)
