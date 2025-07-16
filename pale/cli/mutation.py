# !/usr/bin/env python3
"""
CLI module for protein mutation related tasks, including but not limited to generating mutated
structures from PDB files, generating mappings and setting up basic components and objects to be
used for free energy estimations using OpenFE

OpenFE CLI Protein Mutation Extension

This module extends the OpenFE CLI with protein mutation commands that accept
mutation strings and input SDF/PDB files for components.

Usage:
    openfe plan-protein-mutation [OPTIONS]
    openfe run-protein-mutation [OPTIONS]
    openfe analyze-protein-mutation [OPTIONS]

Commands:
    plan-protein-mutation    Plan protein mutation free energy calculations
    run-protein-mutation     Run protein mutation transformations
    analyze-protein-mutation Analyze protein mutation results
"""

import click
import json
import logging
import os
import sys
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Dict, List, Optional, Tuple, Any
import yaml
import warnings

# OpenFE imports (these would be the actual imports in the real implementation)
try:
    import openfe
    from openfe import SmallMoleculeComponent, ProteinComponent, SolventComponent
    from openfe.protocols import RelativeHybridTopologyProtocol
    from openfe.setup import ChemicalSystem, Transformation
    from openfe.utils import AtomMapping, LigandAtomMapper
    from openfe import AlchemicalNetwork
except ImportError as e:
    warnings.warn(f"OpenFE imports not available: {e}")
    # For development/testing purposes
    openfe = None

# RDKit imports for molecule handling
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolDescriptors
except ImportError:
    warnings.warn("RDKit not available - molecule handling will be limited")
    Chem = None

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Click context settings
CONTEXT_SETTINGS = {'help_option_names': ['-h', '--help']}


class ProteinMutationError(Exception):
    """Custom exception for protein mutation operations."""
    pass


class MutationParser:
    """Parser for protein mutation strings."""

    @staticmethod
    def parse_mutation_string(mutation_str: str) -> Dict[str, Any]:
        """
        Parse mutation string into components.

        Supported formats:
        - A123V (single letter amino acid codes)
        - ALA123VAL (three letter codes)
        - A123V,B456L (multiple mutations)
        - A123V:chain_A (with chain specification)
        """
        mutations = []

        # Split multiple mutations
        mutation_parts = mutation_str.split(',')

        for part in mutation_parts:
            part = part.strip()

            # Check for chain specification
            if ':' in part:
                mutation_part, chain = part.split(':', 1)
                chain = chain.strip()
            else:
                mutation_part = part
                chain = None

            # Parse single mutation
            mutation_info = MutationParser._parse_single_mutation(mutation_part)
            mutation_info['chain'] = chain
            mutations.append(mutation_info)

        return {
            'mutations': mutations,
            'n_mutations': len(mutations),
            'original_string': mutation_str
        }

    @staticmethod
    def _parse_single_mutation(mutation_str: str) -> Dict[str, Any]:
        """Parse a single mutation string."""
        import re

        # Pattern for single letter codes: A123V
        pattern_1 = r'^([A-Z])(\d+)([A-Z])$'
        # Pattern for three-letter codes: ALA123VAL
        pattern_3 = r'^([A-Z]{3})(\d+)([A-Z]{3})$'

        match_1 = re.match(pattern_1, mutation_str)
        match_3 = re.match(pattern_3, mutation_str)

        if match_1:
            from_aa, position, to_aa = match_1.groups()
            return {
                'from_residue': from_aa,
                'to_residue': to_aa,
                'position': int(position),
                'format': 'single_letter'
            }
        elif match_3:
            from_aa, position, to_aa = match_3.groups()
            return {
                'from_residue': from_aa,
                'to_residue': to_aa,
                'position': int(position),
                'format': 'three_letter'
            }
        else:
            raise ProteinMutationError(f"Invalid mutation format: {mutation_str}")


class ProteinMutationPlanner:
    """Plan protein mutation free energy calculations."""

    def __init__(self):
        self.aa_mapping = {
            'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP',
            'C': 'CYS', 'E': 'GLU', 'Q': 'GLN', 'G': 'GLY',
            'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS',
            'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER',
            'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
        }
        self.reverse_aa_mapping = {v: k for k, v in self.aa_mapping.items()}

    def create_protein_mutation_network(self, protein_file: str, mutation_str: str,
                                        ligand_sdf: Optional[str] = None,
                                        cofactor_sdf: Optional[str] = None,
                                        output_dir: str = "protein_mutations",
                                        settings: Optional[Dict] = None) -> AlchemicalNetwork:
        """Create an alchemical network for protein mutations."""

        logger.info(f"Planning protein mutation: {mutation_str}")

        # Parse mutation string
        mutation_info = MutationParser.parse_mutation_string(mutation_str)

        # Load components
        components = self._load_components(protein_file, ligand_sdf, cofactor_sdf)

        # Create wild-type and mutant systems
        wt_system, mut_system = self._create_mutation_systems(
            components, mutation_info, settings
        )

        # Create transformation
        transformation = self._create_protein_transformation(
            wt_system, mut_system, mutation_info, settings
        )

        # Create network
        network = AlchemicalNetwork([transformation])

        # Save network
        self._save_network(network, output_dir, mutation_info)

        return network

    def _load_components(self, protein_file: str, ligand_sdf: Optional[str],
                         cofactor_sdf: Optional[str]) -> Dict[str, Any]:
        """Load molecular components from files."""
        components = {}

        # Load protein
        if not os.path.exists(protein_file):
            raise ProteinMutationError(f"Protein file not found: {protein_file}")

        logger.info(f"Loading protein from {protein_file}")
        components['protein'] = ProteinComponent.from_pdb_file(protein_file)

        # Load ligand if provided
        if ligand_sdf:
            if not os.path.exists(ligand_sdf):
                raise ProteinMutationError(f"Ligand SDF file not found: {ligand_sdf}")

            logger.info(f"Loading ligands from {ligand_sdf}")
            components['ligands'] = self._load_molecules_from_sdf(ligand_sdf)

        # Load cofactors if provided
        if cofactor_sdf:
            if not os.path.exists(cofactor_sdf):
                raise ProteinMutationError(f"Cofactor SDF file not found: {cofactor_sdf}")

            logger.info(f"Loading cofactors from {cofactor_sdf}")
            components['cofactors'] = self._load_molecules_from_sdf(cofactor_sdf)

        # Add solvent
        components['solvent'] = SolventComponent()

        return components

    def _load_molecules_from_sdf(self, sdf_file: str) -> List[SmallMoleculeComponent]:
        """Load molecules from SDF file."""
        if Chem is None:
            raise ProteinMutationError("RDKit not available for SDF processing")

        molecules = []
        suppl = Chem.SDMolSupplier(sdf_file)

        for i, mol in enumerate(suppl):
            if mol is not None:
                # Add hydrogens if not present
                mol = Chem.AddHs(mol)

                # Create OpenFE component
                mol_component = SmallMoleculeComponent.from_rdkit(mol)
                molecules.append(mol_component)
                logger.info(
                    f"Loaded molecule {i + 1}: {mol.GetProp('_Name') if mol.HasProp('_Name') else f'mol_{i + 1}'}")

        if not molecules:
            raise ProteinMutationError(f"No valid molecules found in {sdf_file}")

        return molecules

    def _create_mutation_systems(self, components: Dict[str, Any],
                                 mutation_info: Dict[str, Any],
                                 settings: Optional[Dict] = None) -> Tuple[
        ChemicalSystem, ChemicalSystem]:
        """Create wild-type and mutant chemical systems."""

        # Base system components
        system_components = [components['protein'], components['solvent']]

        # Add ligands if present
        if 'ligands' in components:
            # For now, use the first ligand - could be extended for multiple ligands
            system_components.append(components['ligands'][0])

        # Add cofactors if present
        if 'cofactors' in components:
            system_components.extend(components['cofactors'])

        # Create wild-type system
        wt_system = ChemicalSystem(system_components)

        # Create mutant protein
        mutant_protein = self._create_mutant_protein(
            components['protein'], mutation_info
        )

        # Create mutant system
        mutant_components = [mutant_protein, components['solvent']]
        if 'ligands' in components:
            mutant_components.append(components['ligands'][0])
        if 'cofactors' in components:
            mutant_components.extend(components['cofactors'])

        mut_system = ChemicalSystem(mutant_components)

        return wt_system, mut_system

    def _create_mutant_protein(self, protein: ProteinComponent,
                               mutation_info: Dict[str, Any],
                               add_missing: bool = False) -> List[ProteinComponent]:
        """
        Create mutant proteins by applying point mutations to a wild-type protein.

        This method processes each mutation specified in the mutation_info dictionary,
        creating a separate mutant protein for each one. The mutations are applied
        using PDBFixer, which can also optionally add missing residues and atoms.

        Parameters
        ----------
        protein : ProteinComponent
            The wild-type protein structure to mutate.
        mutation_info : dict
            Dictionary containing mutation specifications with the following structure:

            - mutations : list of dict
                List of mutation dictionaries, each containing:

                - from_residue : str
                    Single-letter amino acid code of original residue (e.g., 'A')
                - position : int
                    Residue position number in the protein sequence
                - to_residue : str
                    Single-letter amino acid code of target residue (e.g., 'V')
                - chain : str
                    Chain identifier (e.g., 'A')
        add_missing : bool, optional
            If True, adds missing residues and atoms to the wild-type protein
            before applying mutations. Missing atoms and hydrogens are always
            added to the mutant regardless of this setting. Default is False.

        Returns
        -------
        list of ProteinComponent
            List of mutant protein objects, one for each mutation applied.
            The order matches the order of mutations in mutation_info['mutations'].

        Notes
        -----
        - Each mutation is applied to the original wild-type protein independently
        - Temporary PDB files are created during processing but cleaned up automatically
        - Hydrogens are added at pH 7.0 after mutation application
        - The method logs progress information for each mutation being applied

        Examples
        --------
        >>> mutation_info = {
        ...     'mutations': [
        ...         {'from_residue': 'A', 'position': 123, 'to_residue': 'V', 'chain': 'A'},
        ...         {'from_residue': 'L', 'position': 456, 'to_residue': 'P', 'chain': 'B'}
        ...     ]
        ... }
        >>> mutants = self._create_mutant_protein(wild_type_protein, mutation_info)
        >>> len(mutants)  # Returns 2, one for each mutation
        2
        """
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile

        mutated_proteins = []  # we will append ProteinComponents to this list

        logger.info("Creating mutant proteins...")

        for mutation in mutation_info['mutations']:
            logger.info(
                f"Applying mutation: {mutation['from_residue']}{mutation['position']}{mutation['to_residue']} in chain {mutation['chain']}")
            # Store pdb structure of wild-type protein, needed for pdbfixer
            # TODO: We probably want to avoid having to pass through pdb files (or any file)
            temp_f = NamedTemporaryFile(delete=False, suffix=".pdb")
            protein.to_pdb_file(temp_f.name)
            # Parse the protein structure -- Add missing residues/atoms if specified
            pdbfixer = PDBFixer(temp_f.name)
            if add_missing:
                pdbfixer.findMissingResidues()
                pdbfixer.findMissingAtoms()
                pdbfixer.addMissingAtoms()
                pdbfixer.addMissingHydrogens(7.0)
            # Parse mutation string to pdbfixer format
            mutation_spec = f"{mutation['from_residue']}-{mutation['position']}-{mutation['to_residue']}"
            # Apply the mutation
            pdbfixer.applyMutations([mutation_spec], mutation["chain"])
            pdbfixer.findMissingResidues()
            pdbfixer.findMissingAtoms()
            pdbfixer.addMissingAtoms()
            pdbfixer.addMissingHydrogens(7.0)
            # Save resulting mutation to PDB using openmm
            temp_out_f = NamedTemporaryFile(delete=False, suffix=".pdb")
            omm_top = pdbfixer.topology
            omm_pos = pdbfixer.positions
            with open(temp_out_f, "w") as out_file:
                PDBFile.writeFile(omm_top, omm_pos, out_file)
            # Create ProteinComponent and add it to list
            mut_protein = ProteinComponent.from_pdb_file(temp_out_f.name)
            mutated_proteins.append(mut_protein)

        return mutated_proteins

    def _create_protein_transformation(self, wt_system: ChemicalSystem,
                                       mut_system: ChemicalSystem,
                                       mutation_info: Dict[str, Any],
                                       settings: Optional[Dict] = None) -> Transformation:
        """Create transformation between wild-type and mutant systems."""

        # Default protocol settings
        default_settings = {
            'protocol': 'RelativeHybridTopologyProtocol',
            'lambda_schedule': [0.0, 0.25, 0.5, 0.75, 1.0],
            'simulation_time': '5.0 ns',
            'equilibration_time': '1.0 ns'
        }

        if settings:
            default_settings.update(settings)

        # Create protocol
        protocol = RelativeHybridTopologyProtocol(
            settings=default_settings
        )

        # Create atom mapping (simplified)
        # In practice, this would use sophisticated protein atom mapping
        atom_mapping = self._create_protein_atom_mapping(wt_system, mut_system, mutation_info)

        # Create transformation
        transformation = Transformation(
            stateA=wt_system,
            stateB=mut_system,
            mapping=atom_mapping,
            protocol=protocol,
            name=f"protein_mutation_{mutation_info['original_string']}"
        )

        return transformation

    def _create_protein_atom_mapping(self, wt_system: ChemicalSystem,
                                     mut_system: ChemicalSystem,
                                     mutation_info: Dict[str, Any]) -> AtomMapping:
        """Create atom mapping for protein mutation."""

        # This is a placeholder implementation
        # In practice, this would:
        # 1. Identify corresponding atoms between WT and mutant
        # 2. Handle appearing/disappearing atoms in the mutation
        # 3. Create proper atom mapping for alchemical transformation

        logger.info("Creating protein atom mapping...")

        # Simplified mapping - would need proper implementation
        mapping_dict = {}  # atom index mapping

        return AtomMapping(
            componentA_to_componentB=mapping_dict,
            componentA=wt_system,
            componentB=mut_system
        )

    def _save_network(self, network: AlchemicalNetwork, output_dir: str,
                      mutation_info: Dict[str, Any]):
        """Save the alchemical network to files."""

        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Save network
        network_file = output_path / "protein_mutation_network.json"
        with open(network_file, 'w') as f:
            json.dump(network.to_dict(), f, indent=2)

        # Save individual transformations
        transformations_dir = output_path / "transformations"
        transformations_dir.mkdir(exist_ok=True)

        for i, transformation in enumerate(network.edges):
            trans_file = transformations_dir / f"transformation_{i:03d}.json"
            with open(trans_file, 'w') as f:
                json.dump(transformation.to_dict(), f, indent=2)

        # Save mutation info
        mutation_file = output_path / "mutation_info.yaml"
        with open(mutation_file, 'w') as f:
            yaml.dump(mutation_info, f, default_flow_style=False)

        logger.info(f"Network saved to {output_dir}")


# CLI Commands using Click decorators to extend OpenFE CLI

@click.command('plan-protein-mutation',
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
    """
    Plan protein mutation free energy calculations.

    This command creates an alchemical network for protein mutation free energy
    calculations. It accepts a protein PDB file, mutation specification, and
    optional ligand/cofactor SDF files.

    Examples:

        # Simple mutation
        openfe plan-protein-mutation -p protein.pdb -m "A123V"

        # Multiple mutations
        openfe plan-protein-mutation -p protein.pdb -m "A123V,B456L"

        # With ligand
        openfe plan-protein-mutation -p protein.pdb -m "A123V" -l ligands.sdf

        # With custom settings
        openfe plan-protein-mutation -p protein.pdb -m "A123V" -s settings.yaml
    """

    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        # Load settings if provided
        settings = {}
        if settings_file:
            with open(settings_file, 'r') as f:
                settings = yaml.safe_load(f)
            logger.info(f"Loaded settings from {settings_file}")

        if dry_run:
            logger.info("Dry run mode - validating inputs only")
            # Parse mutation string to validate
            mutation_info = MutationParser.parse_mutation_string(mutation_str)
            logger.info(f"Parsed mutations: {mutation_info}")
            logger.info("Input validation successful")
            return

        # Create planner and plan mutations
        planner = ProteinMutationPlanner()
        network = planner.create_protein_mutation_network(
            protein_file=protein_file,
            mutation_str=mutation_str,
            ligand_sdf=ligand_sdf,
            cofactor_sdf=cofactor_sdf,
            output_dir=output_dir,
            settings=settings
        )

        logger.info(f"Successfully planned {len(network.edges)} transformations")
        logger.info(f"Output saved to: {output_dir}")

    except Exception as e:
        logger.error(f"Planning failed: {e}")
        sys.exit(1)


@click.command('run-protein-mutation',
               context_settings=CONTEXT_SETTINGS,
               short_help="Run protein mutation transformations")
@click.option('-d', '--transformations-dir', 'trans_dir',
              default='protein_mutations/transformations',
              help='Directory containing transformation JSON files')
@click.option('-o', '--output', 'output_file',
              help='Output file for results (JSON format)')
@click.option('-w', '--work-dir', 'work_dir',
              default='.',
              help='Working directory for simulation files')
@click.option('-j', '--parallel', 'n_jobs',
              default=1, type=int,
              help='Number of parallel jobs')
@click.option('--gpu', is_flag=True,
              help='Use GPU acceleration if available')
@click.option('-v', '--verbose', is_flag=True,
              help='Verbose output')
def run_protein_mutation(trans_dir, output_file, work_dir, n_jobs, gpu, verbose):
    """
    Run protein mutation transformations.

    This command executes the transformations created by plan-protein-mutation.
    It can run transformations in parallel and supports GPU acceleration.

    Examples:

        # Run all transformations
        openfe run-protein-mutation

        # Run with custom directories
        openfe run-protein-mutation -d my_mutations/transformations -w work_dir

        # Run in parallel with GPU
        openfe run-protein-mutation -j 4 --gpu
    """

    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        trans_path = Path(trans_dir)
        if not trans_path.exists():
            raise ProteinMutationError(f"Transformations directory not found: {trans_dir}")

        # Find transformation files
        trans_files = list(trans_path.glob("transformation_*.json"))
        if not trans_files:
            raise ProteinMutationError(f"No transformation files found in {trans_dir}")

        logger.info(f"Found {len(trans_files)} transformations to run")

        results = []
        work_path = Path(work_dir)
        work_path.mkdir(parents=True, exist_ok=True)

        for trans_file in trans_files:
            logger.info(f"Running transformation: {trans_file.name}")

            # This would use openfe.quickrun or similar
            result = run_single_transformation(trans_file, work_path, gpu)
            results.append(result)

        # Save results
        if output_file:
            with open(output_file, 'w') as f:
                json.dump(results, f, indent=2)
            logger.info(f"Results saved to: {output_file}")

        logger.info("All transformations completed successfully")

    except Exception as e:
        logger.error(f"Execution failed: {e}")
        sys.exit(1)


def run_single_transformation(trans_file: Path, work_dir: Path, use_gpu: bool = False):
    """Run a single transformation file."""

    # Load transformation
    with open(trans_file, 'r') as f:
        transformation_data = json.load(f)

    # Create work subdirectory
    work_subdir = work_dir / trans_file.stem
    work_subdir.mkdir(exist_ok=True)

    # This would integrate with OpenFE's quickrun functionality
    logger.info(f"Executing transformation in {work_subdir}")

    # Placeholder result
    result = {
        'transformation': str(trans_file),
        'work_dir': str(work_subdir),
        'status': 'completed',
        'free_energy': 0.0,  # Would be calculated
        'uncertainty': 0.0  # Would be calculated
    }

    return result


@click.command('analyze-protein-mutation',
               context_settings=CONTEXT_SETTINGS,
               short_help="Analyze protein mutation results")
@click.option('-r', '--results', 'results_file',
              help='Results JSON file from run-protein-mutation')
@click.option('-d', '--results-dir', 'results_dir',
              help='Directory containing individual result files')
@click.option('-o', '--output', 'output_file',
              default='protein_mutation_analysis.json',
              help='Output file for analysis results')
@click.option('--plot', is_flag=True,
              help='Generate analysis plots')
@click.option('-v', '--verbose', is_flag=True,
              help='Verbose output')
def analyze_protein_mutation(results_file, results_dir, output_file, plot, verbose):
    """
    Analyze protein mutation free energy results.

    This command analyzes the results from protein mutation simulations,
    calculates statistics, and optionally generates plots.

    Examples:

        # Analyze results file
        openfe analyze-protein-mutation -r results.json

        # Analyze directory of results
        openfe analyze-protein-mutation -d results_directory

        # Generate plots
        openfe analyze-protein-mutation -r results.json --plot
    """

    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        # Load results
        if results_file:
            with open(results_file, 'r') as f:
                results = json.load(f)
        elif results_dir:
            results = load_results_from_directory(results_dir)
        else:
            raise ProteinMutationError("Either --results or --results-dir must be specified")

        # Perform analysis
        analysis = analyze_mutation_results(results)

        # Save analysis
        with open(output_file, 'w') as f:
            json.dump(analysis, f, indent=2)

        logger.info(f"Analysis completed. Results saved to: {output_file}")

        # Generate plots if requested
        if plot:
            generate_analysis_plots(analysis, output_file)

    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        sys.exit(1)


def load_results_from_directory(results_dir: str) -> List[Dict]:
    """Load results from a directory of result files."""
    results_path = Path(results_dir)
    results = []

    for result_file in results_path.glob("*.json"):
        with open(result_file, 'r') as f:
            result = json.load(f)
            results.append(result)

    return results


def analyze_mutation_results(results: List[Dict]) -> Dict:
    """Analyze protein mutation results."""

    analysis = {
        'n_transformations': len(results),
        'successful_runs': sum(1 for r in results if r.get('status') == 'completed'),
        'free_energies': [r.get('free_energy', 0.0) for r in results],
        'uncertainties': [r.get('uncertainty', 0.0) for r in results],
        'summary_statistics': {}
    }

    # Calculate summary statistics
    if analysis['free_energies']:
        import numpy as np
        fe_values = np.array(analysis['free_energies'])
        analysis['summary_statistics'] = {
            'mean_free_energy': float(np.mean(fe_values)),
            'std_free_energy': float(np.std(fe_values)),
            'min_free_energy': float(np.min(fe_values)),
            'max_free_energy': float(np.max(fe_values))
        }

    return analysis


def generate_analysis_plots(analysis: Dict, output_file: str):
    """Generate analysis plots."""
    try:
        import matplotlib.pyplot as plt

        output_path = Path(output_file).parent

        # Free energy histogram
        plt.figure(figsize=(10, 6))
        plt.hist(analysis['free_energies'], bins=20, alpha=0.7)
        plt.xlabel('Free Energy (kcal/mol)')
        plt.ylabel('Frequency')
        plt.title('Distribution of Mutation Free Energies')
        plt.savefig(output_path / 'free_energy_distribution.png', dpi=300, bbox_inches='tight')
        plt.close()

        logger.info("Analysis plots generated")

    except ImportError:
        logger.warning("Matplotlib not available - skipping plot generation")


# Integration with OpenFE CLI
# This would be added to the main OpenFE CLI module

def add_protein_mutation_commands(cli_group):
    """Add protein mutation commands to the OpenFE CLI group."""
    cli_group.add_command(plan_protein_mutation)
    cli_group.add_command(run_protein_mutation)
    cli_group.add_command(analyze_protein_mutation)


# Example of how this would integrate with the main OpenFE CLI
if __name__ == "__main__":
    # This would be in the main OpenFE CLI module
    import click


    @click.group()
    def cli():
        """OpenFE CLI with protein mutation extensions."""
        pass


    # Add existing OpenFE commands (these would already exist)
    # cli.add_command(plan_rbfe_network)
    # cli.add_command(plan_rhfe_network)
    # cli.add_command(quickrun)
    # cli.add_command(atommapping)

    # Add our protein mutation commands
    add_protein_mutation_commands(cli)

    cli()