"""
Config file for tests. Place to store some base fixture, pytest options or modifiers.
"""

from importlib.resources import files
import pytest


@pytest.fixture(scope="session")
def test_data_path():
    data_path = files("pale.tests.data")
    return data_path


@pytest.fixture(scope="session")
def dipeptide_test_path(test_data_path):
    file_path = test_data_path.joinpath("ALA_capped.pdb")
    return file_path


@pytest.fixture(scope="session")
def dipeptide_mutate_spec(dipeptide_test_path):
    """
    Generate a mutation specification for the central residue of a dipeptide.

    The function reads a PDB file of a dipeptide (two peptide bonds, three amino acids: X–Y–Z),
    autodetects the chain ID, and returns a mutation specification string split into
    the residue mutation part and the chain ID.

    Parameters
    ----------
    dipeptide_test_path : pathlib.Path
        Path to the PDB file containing the dipeptide structure.

    Returns
    -------
    list of str
        Two-element list containing:
        1. Mutation specification in the format ``<current_aa>-<residue_index>-<mutate_aa>``.
        2. Chain ID of the residue to mutate.

    Notes
    -----
    - The mutation is always from the second residue in the sequence to proline ("PRO").
    - Residue numbering is assumed to start at 1.
    """
    from openmm.app import PDBFile

    pdb_file = PDBFile(dipeptide_test_path)
    omm_top = pdb_file.topology
    residue = list(omm_top.residues())[1]  # Extract 2nd residue
    chain = residue.chain.id  # Extract chain id/name
    current_aa = residue.name
    mutate_aa = "PRO"
    mut_spec = f"{current_aa}-2-{mutate_aa}:{chain}"
    return mut_spec.split(":")
