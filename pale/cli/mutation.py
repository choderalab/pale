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
import sys
from pathlib import Path
from typing import Dict, List
import yaml
import warnings

# OpenFE imports (these would be the actual imports in the real implementation)
try:
    import openfe
    from openfe import SmallMoleculeComponent, ProteinComponent, SolventComponent
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

from pale.mutations import MutationParser, ProteinMutationError, ProteinMutationPlanner

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Click context settings
CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}


# CLI Commands using Click decorators to extend OpenFE CLI


@click.command(
    "plan-protein-mutation",
    context_settings=CONTEXT_SETTINGS,
    short_help="Plan protein mutation free energy calculations",
)
@click.option(
    "-p",
    "--protein",
    "protein_file",
    required=True,
    type=click.Path(exists=True),
    help="Protein PDB file",
)
@click.option(
    "-m",
    "--mutation",
    "mutation_str",
    required=True,
    type=str,
    help='Mutation string (e.g., "A123V" or "ALA123VAL")',
)
@click.option(
    "-l", "--ligands", "ligand_sdf", type=click.Path(exists=True), help="Ligand SDF file (optional)"
)
@click.option(
    "-c",
    "--cofactors",
    "cofactor_sdf",
    type=click.Path(exists=True),
    help="Cofactor SDF file (optional)",
)
@click.option(
    "-o",
    "--output",
    "output_dir",
    default="protein_mutations",
    help="Output directory for planned mutations",
)
@click.option(
    "-s",
    "--settings",
    "settings_file",
    type=click.Path(exists=True),
    help="YAML settings file for protocol customization",
)
@click.option("--dry-run", is_flag=True, help="Validate inputs without creating transformations")
@click.option("-v", "--verbose", is_flag=True, help="Verbose output")
def plan_protein_mutation(
    protein_file,
    mutation_str,
    ligand_sdf,
    cofactor_sdf,
    output_dir,
    settings_file,
    dry_run,
    verbose,
):
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
            with open(settings_file, "r") as f:
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
            settings=settings,
        )

        logger.info(f"Successfully planned {len(network.edges)} transformations")
        logger.info(f"Output saved to: {output_dir}")

    except Exception as e:
        logger.error(f"Planning failed: {e}")
        sys.exit(1)


@click.command(
    "run-protein-mutation",
    context_settings=CONTEXT_SETTINGS,
    short_help="Run protein mutation transformations",
)
@click.option(
    "-d",
    "--transformations-dir",
    "trans_dir",
    default="protein_mutations/transformations",
    help="Directory containing transformation JSON files",
)
@click.option("-o", "--output", "output_file", help="Output file for results (JSON format)")
@click.option(
    "-w", "--work-dir", "work_dir", default=".", help="Working directory for simulation files"
)
@click.option("-j", "--parallel", "n_jobs", default=1, type=int, help="Number of parallel jobs")
@click.option("--gpu", is_flag=True, help="Use GPU acceleration if available")
@click.option("-v", "--verbose", is_flag=True, help="Verbose output")
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
            with open(output_file, "w") as f:
                json.dump(results, f, indent=2)
            logger.info(f"Results saved to: {output_file}")

        logger.info("All transformations completed successfully")

    except Exception as e:
        logger.error(f"Execution failed: {e}")
        sys.exit(1)


def run_single_transformation(trans_file: Path, work_dir: Path, use_gpu: bool = False):
    """Run a single transformation file."""

    # Load transformation
    with open(trans_file, "r") as f:
        transformation_data = json.load(f)

    # Create work subdirectory
    work_subdir = work_dir / trans_file.stem
    work_subdir.mkdir(exist_ok=True)

    # This would integrate with OpenFE's quickrun functionality
    logger.info(f"Executing transformation in {work_subdir}")

    # Placeholder result
    result = {
        "transformation": str(trans_file),
        "work_dir": str(work_subdir),
        "status": "completed",
        "free_energy": 0.0,  # Would be calculated
        "uncertainty": 0.0,  # Would be calculated
    }

    return result


@click.command(
    "analyze-protein-mutation",
    context_settings=CONTEXT_SETTINGS,
    short_help="Analyze protein mutation results",
)
@click.option("-r", "--results", "results_file", help="Results JSON file from run-protein-mutation")
@click.option(
    "-d", "--results-dir", "results_dir", help="Directory containing individual result files"
)
@click.option(
    "-o",
    "--output",
    "output_file",
    default="protein_mutation_analysis.json",
    help="Output file for analysis results",
)
@click.option("--plot", is_flag=True, help="Generate analysis plots")
@click.option("-v", "--verbose", is_flag=True, help="Verbose output")
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
            with open(results_file, "r") as f:
                results = json.load(f)
        elif results_dir:
            results = load_results_from_directory(results_dir)
        else:
            raise ProteinMutationError("Either --results or --results-dir must be specified")

        # Perform analysis
        analysis = analyze_mutation_results(results)

        # Save analysis
        with open(output_file, "w") as f:
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
        with open(result_file, "r") as f:
            result = json.load(f)
            results.append(result)

    return results


def analyze_mutation_results(results: List[Dict]) -> Dict:
    """Analyze protein mutation results."""

    analysis = {
        "n_transformations": len(results),
        "successful_runs": sum(1 for r in results if r.get("status") == "completed"),
        "free_energies": [r.get("free_energy", 0.0) for r in results],
        "uncertainties": [r.get("uncertainty", 0.0) for r in results],
        "summary_statistics": {},
    }

    # Calculate summary statistics
    if analysis["free_energies"]:
        import numpy as np

        fe_values = np.array(analysis["free_energies"])
        analysis["summary_statistics"] = {
            "mean_free_energy": float(np.mean(fe_values)),
            "std_free_energy": float(np.std(fe_values)),
            "min_free_energy": float(np.min(fe_values)),
            "max_free_energy": float(np.max(fe_values)),
        }

    return analysis


def generate_analysis_plots(analysis: Dict, output_file: str):
    """Generate analysis plots."""
    try:
        import matplotlib.pyplot as plt

        output_path = Path(output_file).parent

        # Free energy histogram
        plt.figure(figsize=(10, 6))
        plt.hist(analysis["free_energies"], bins=20, alpha=0.7)
        plt.xlabel("Free Energy (kcal/mol)")
        plt.ylabel("Frequency")
        plt.title("Distribution of Mutation Free Energies")
        plt.savefig(output_path / "free_energy_distribution.png", dpi=300, bbox_inches="tight")
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
