"""
Test for the CLI in pale

Aimed towards running hopefully all the CLI commands and options and check that they return
expected values.
"""

from click.testing import CliRunner
from pale.cli import mutate


def test_mutate_from_pdb(dipeptide_test_path, dipeptide_mutate_spec):
    """
    Test cli to generate protein mutation PDB
    """
    runner = CliRunner()
    result = runner.invoke(
        mutate, args=f"--input-pdb {dipeptide_test_path} --mutation-spec {dipeptide_mutate_spec}"
    )
    assert result.exit_code == 0, (
        f"Command returned unsuccessfully with return code {result.exit_code}."
    )


def test_failed_mutate_from_pdb():
    """
    Test CLI failing for generate protein mutation.
    """
    runner = CliRunner()
    result = runner.invoke(mutate, args=f"--input-pdb {test_pdb} --mutation-spec {mut_spec}")
    assert result.exit_code != 0, "Exit code was expected to be different than zero for this test."
