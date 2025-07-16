"""
Module with tools to run specific units from a `gufe` protocol DAG.
"""

import os
import click

@click.command('run-unit',
               context_settings=CONTEXT_SETTINGS,
               short_help="Plan protein mutation free energy calculations")
@click.option('-p', '--protein', 'protein_file',
              required=True, type=click.Path(exists=True),
              help='Protein PDB file')
@click.option('-m', '--mutation', 'mutation_str',
              required=True, type=str,
              help='Mutation string (e.g., "A123V" or "ALA123VAL")')
@click.option('-l', '--ligands', 'ligand_sdf',
              type=click.Path(exists=True),
              help='Ligand SDF file (optional)')
@click.option('-c', '--cofactors', 'cofactor_sdf',
              type=click.Path(exists=True),
              help='Cofactor SDF file (optional)')
@click.option('-o', '--output', 'output_dir',
              default='protein_mutations',
              help='Output directory for planned mutations')
@click.option('-s', '--settings', 'settings_file',
              type=click.Path(exists=True),
              help='YAML settings file for protocol customization')
@click.option('--dry-run', is_flag=True,
              help='Validate inputs without creating transformations')
@click.option('-v', '--verbose', is_flag=True,
              help='Verbose output')
def plan_protein_mutation(protein_file, mutation_str, ligand_sdf, cofactor_sdf,
                          output_dir, settings_file, dry_run, verbose):
    pass


def run_cycle_unit(system_path: str | os.PathLike, state_path: str | os.PathLike,
                   integrator_path: str | os.PathLike, htf_path: str | os.PathLike):
    from feflow.utils.data import deserialize
    import numpy as np

    # Deserialize OpenMM objects
    system = deserialize(system_path)
    state = deserialize(state_path)
    integrator = deserialize(integrator_path)

    # Extract HTF information from pickled object
    hybrid_top_factory = np.load(htf_path)
    initial_atom_indices = hybrid_top_factory.initial_atom_indices
    final_atom_indices = hybrid_top_factory.final_atom_indices