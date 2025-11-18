#!/usr/bin/env python
# coding: utf-8

# Modified from Iván Pulido + OpenFe docs

import argparse
import os
import pathlib
from pathlib import Path
import re
import pandas as pd
import pdbfixer
from openmm.app import PDBFile


##########################################################################################################

# Read args
parser = argparse.ArgumentParser(
    description="Run DAG protocol for protein mutations using PALE."
)
parser.add_argument(
    "-p",
    dest="protein",
    type=str,
    help="Path to input protein PDB file"
)

parser.add_argument(
    "-l",
    dest="ligand",
    type=str,
    default="",
    help="Optional path to input ligand SDF file"
)

parser.add_argument(
    "-m",
    dest="mutation",
    type=str,
    help="Mutation using one letter code for residue (i.e. L1A)"
)

parser.add_argument(
    "--leg",
    dest="leg",
    type=str,
    help="Thermodynamic cycle leg to be added as prefix for output files"
)

parser.add_argument(
    "--mutant_chain",
    dest="mutant_chain",
    type=str,
    default="A",
    help="Optional ID of chain with mutant residue in PDB file"
)

parser.add_argument(
    "--num_cycles",
    dest="num_cycles",
    type=int,
    help="Number of cycles to be run for DAG"
)
args = parser.parse_args()

##########################################################################################################
# Util code

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
    "GLN": "Q"
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
    pattern = r'([A-Z]{1,3})(\d+)([A-Z]{1,3})'
    match = re.search(pattern, mutation_string)
    if match:
        initial_aa, res_number, final_aa = match.groups()
        return initial_aa, res_number, final_aa
    raise ValueError(f"Invalid mutation specification: {mutation_spec}")

##########################################################################################################
##########################################################################################################

# initiate pdb_fixer for each mutation
pdb_fixer = pdbfixer.PDBFixer(filename=args.protein)
pdb_fixer.findMissingResidues()
omm_top = pdb_fixer.topology # store topology

# Get mutation 
initial_aa, res_num, final_aa = parse_mutation_spec(args.mutation)
mutation_str = "-".join([aa_one_to_three_code[initial_aa], res_num, aa_one_to_three_code[final_aa]])

# Apply mutation
pdb_fixer.applyMutations(mutations=[mutation_str], chain_id=args.mutant_chain)
pdb_fixer.findMissingResidues()
pdb_fixer.findMissingAtoms()
pdb_fixer.addMissingAtoms()
pdb_fixer.addMissingHydrogens()
omm_topology = pdb_fixer.topology
omm_positions = pdb_fixer.positions
mutant_out_dir = f"./abl/{args.leg}/mutant_structures"
os.makedirs(mutant_out_dir, exist_ok=True)
with open(f"{mutant_out_dir}/{args.leg}_{mutation_str}.pdb", "w") as out_file:
    PDBFile.writeFile(omm_topology, omm_positions, out_file)


from kartograf import KartografAtomMapper
from gufe import ProteinComponent, SmallMoleculeComponent, SolventComponent, ChemicalSystem, LigandAtomMapping

# Generating mappings for all relevant mutations
reference_component = ProteinComponent.from_pdb_file(args.protein)
atom_mapper = KartografAtomMapper()

print(f"Generating {mutation_str} mutation mapping.")
# Read mutant pdb component
mutant_comp = ProteinComponent.from_pdb_file(f"{mutant_out_dir}/{args.leg}_{mutation_str}.pdb")
# Generate mappings
mapping = next(atom_mapper.suggest_mappings(reference_component, mutant_comp))
mapping_dir = f"./abl/{args.leg}/mappings"
# Serializing mappings
os.makedirs(mapping_dir, exist_ok=True)
with open(f"{mapping_dir}/{args.leg}_{mutation_str}.json", "w") as out_file:
    mapping.to_json(out_file)

from feflow.protocols import NonEquilibriumCyclingProtocol
from gufe.protocols import execute_DAG

# Specify settings for the protocol
settings = NonEquilibriumCyclingProtocol.default_settings()
settings.integrator_settings.equilibrium_steps = 375000  # Short for debugging
settings.integrator_settings.nonequilibrium_steps = 375000
settings.num_cycles = args.num_cycles
settings.alchemical_settings.explicit_charge_correction = True
settings.partial_charge_settings.number_of_conformers = 1
# Using espaloma for ATP partial charge assignment
settings.partial_charge_settings.partial_charge_method = "espaloma"
settings.partial_charge_settings.off_toolkit_backend = "rdkit"


# Create systems and protocol dags
initial_protein_comp = ProteinComponent.from_pdb_file(args.protein)
solvent_comp = SolventComponent()
protocol = NonEquilibriumCyclingProtocol(settings=settings)


final_protein_comp = ProteinComponent.from_pdb_file(f"{mutant_out_dir}/{args.leg}_{mutation_str}.pdb") # mutant protein
if args.ligand != "": # if ligand provided
    print(f"Initializing system with protein, ligand, and solvent.")
    ligand_comp = SmallMoleculeComponent.from_sdf_file(args.ligand)
    initial_state = ChemicalSystem(components={"protein": initial_protein_comp, "ligand": ligand_comp, "solvent": solvent_comp})
    end_state = ChemicalSystem(components={"protein": final_protein_comp, "ligand": ligand_comp, "solvent": solvent_comp})
else: # if ligand not provided
    print(f"Initializing system with protein and solvent.")
    initial_state = ChemicalSystem(components={"protein": initial_protein_comp, "solvent": solvent_comp})
    end_state = ChemicalSystem(components={"protein": final_protein_comp, "solvent": solvent_comp})

mapping = LigandAtomMapping.from_json(f"{mapping_dir}/{args.leg}_{mutation_str}.json")
dag_name = f"NEQCycDAG_{args.leg}_{mutation_str}"
protocol_dag = protocol.create(stateA=initial_state, stateB=end_state, mapping=mapping, 
                               name=dag_name)
protocol_dag_outdir = f"./abl/{args.leg}/protocol_dags"
os.makedirs(protocol_dag_outdir, exist_ok=True)
protocol_dag.to_json(f"{protocol_dag_outdir}/{dag_name}.json")

protocol_dag_path = f"{protocol_dag_outdir}/{dag_name}.json"
protocol_dag_deserialized = NonEquilibriumCyclingProtocol.from_json(protocol_dag_path)


results_path = pathlib.Path(f"./abl/{args.leg}/results_{mutation_str}")
results_path.mkdir(exist_ok=True)
print(f"Executing protocol dag for {mutation_str} mutation.")
protocol_result_dag = execute_DAG(protocol_dag_deserialized, keep_shared=True, shared_basedir=results_path, scratch_basedir=results_path)
protocol_result_dag.to_json(f"./{results_path}/NEqCycDAG_{mutation_str}.json")
print(f"Finished {mutation_str} mutation.")
