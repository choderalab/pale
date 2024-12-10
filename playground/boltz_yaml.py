"""
Script to help predicting protein-ligand complex structures using boltz.

The purpose is to programmatically build a YAML file that boltz can use to
make predictions. It takes data from a PDB and SDF files for protein and
ligands, respectively.

Meant to be used as a tool to help sanity-checking boltz.
"""

import argparse
import yaml
from Bio.Data.IUPACData import protein_letters_3to1
from pdbfixer import PDBFixer
from openff.toolkit import Molecule


def sequence_from_pdb(pdb_path: str):
    """
    Extract the protein sequence from a PDB file. It uses Openmm pdbfixer to extract the
    information.
    """
    fixer = PDBFixer(pdb_path)

    # Get single-letter sequence from pdbfixer object
    # Traverse all chains -- raise warnings if hit with unknown residue code
    sequence = []
    for residue in fixer.topology.residues():
        try:
            sequence.append(protein_letters_3to1[residue.name.title()])
        except KeyError:  # It should only fail with caps
            print(f"Residue {residue} conversion failed. Continuing.")
            continue
    return sequence


def ligands_data_from_sdf(sdf_file: str):
    """
    Generates a dictionary with the ligands data from an sdf file.

    Ligands indices are the keys.
    Ligands smiles are the values.

    TODO: Use ligands names/ids instead of indices. There's an issue with boltz schema.
    """
    # Get smiles for ligands
    ligands = Molecule.from_file(sdf_file)
    smiles_data = {str(mol_index): molecule.to_smiles(explicit_hydrogens=False) for
                   mol_index, molecule in enumerate(ligands)}
    return smiles_data


def generate_boltz_yaml(protein_sequence: list, ligands_data: dict, protein_id: str, output_yaml: str):
    """
    Writes a YAML file that can be used by boltz, from protein sequence and ligands data
    information.

    Note: The resulting YAML files are intended for use with the remote MSA server's capabilities.
    That is, --use_msa_server flag for boltz.
    """
    yaml_dict = {"version": 1, "sequences": []}  # initialize obj to populate
    # Populate protein sequence
    yaml_dict["sequences"].append(
        {"protein": {"id": protein_id, "sequence": "".join(protein_sequence)}})
    # Populate ligands smiles
    for ligand_id, ligand_smiles in ligands_data.items():
        yaml_dict["sequences"].append({"ligand": {"id": ligand_id, "smiles": ligand_smiles}})
    # Dump data in yaml file
    with open(output_yaml, "w") as out_file:
        yaml.dump(yaml_dict, out_file)


def arg_parser():
    parser = argparse.ArgumentParser(
        description="CLI to write YAML files for boltz protein-ligand predictions.")
    parser.add_argument("--protein-file", type=str, help="Path to protein PDB file.")
    parser.add_argument("--ligand-file", type=str, help="Path to liband SDF file.")
    parser.add_argument("--protein-id", type=str, help="Protein name/id.", default="protein")
    parser.add_argument("--output-yaml", type=str, help="Path to yaml output.",
                        default="ligand.yaml")
    return parser


def main():
    parser = arg_parser()
    args = parser.parse_args()
    protein_sequence = sequence_from_pdb(args.protein_file)
    ligands_data = ligands_data_from_sdf(args.ligand_file)
    # Write yaml file
    generate_boltz_yaml(protein_sequence, ligands_data, args.protein_id, args.output_yaml)


if __name__ == "__main__":
    main()
