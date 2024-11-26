"""
Script to predict tyk2 complex structures using boltz.

It takes data from a local copy of the protein ligand benchmark repository
hosted in https://github.com/openforcefield/protein-ligand-benchmark

The purpose is to programatically build a YAML file that boltz can use to
make predictions from.

Meant to be used as a sanity check for boltz.

TODOs:
 * Make basic CLI that accepts paths to pdb/sdf files.
 * Actually generate the YAML file.
"""

import yaml
from Bio.Data.IUPACData import protein_letters_3to1
from pdbfixer import PDBFixer
from openff.toolkit import Molecule

# Read pdb file
pdb_path = "/home/ijpulidos/workdir/repos/protein-ligand-benchmark/data/tyk2/01_protein/crd/protein.pdb"
fixer = PDBfixer(pdb_path)

# Get single-letter sequence from pdbfixer object
sequence_ojb = fixer.sequences[0]  # get first sequence in fixer obj
sequence = []
for residue in sequence_obj.residues:
    try:
        sequence.append(protein_letters_3to1[residue.title()])
    except KeyError:  # It should only fail with caps
        print(f"Residue {residue} conversion failed. Continuing.")
        continue

# Get smiles for ligands
ligands = Molecule.from_file("/home/ijpulidos/workdir/repos/protein-ligand-benchmark/data/tyk2/02_ligands/ligands.sdf")
smiles_data = {molecule.name: molecule.to_smiles(explicit_hydrogens=False) for molecule in ligands}
