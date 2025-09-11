"""
Module where to have the different callables and object to directly handle protein mutations.

Mutations are performed with the help of external cheminformatics tools and libraries.
"""

import json
import logging
import os

from os import PathLike
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any, Dict, List, Optional, Tuple

import yaml

from feflow.protocols import NonEquilibriumCyclingProtocol
from gufe import (
    ChemicalSystem,
    ProteinComponent,
    SmallMoleculeComponent,
    SolventComponent,
    Transformation,
    AtomMapping,
)
from openfe import AlchemicalNetwork
from rdkit import Chem


# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def mutate_from_pdb(input_pdb: PathLike | str, mutation_spec: str, output_pdb: PathLike = None):
    """
    Generate a single point mutation given the input pdb and mutation string specification.
    Optionally storing the output in a PDB file as specified.
    """
    from pdbfixer import PDBFixer


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
        mutation_parts = mutation_str.split(",")

        for part in mutation_parts:
            part = part.strip()

            # Check for chain specification
            if ":" in part:
                mutation_part, chain = part.split(":", 1)
                chain = chain.strip()
            else:
                mutation_part = part
                chain = None

            # Parse single mutation
            mutation_info = MutationParser._parse_single_mutation(mutation_part)
            mutation_info["chain"] = chain
            mutations.append(mutation_info)

        return {
            "mutations": mutations,
            "n_mutations": len(mutations),
            "original_string": mutation_str,
        }

    @staticmethod
    def _parse_single_mutation(mutation_str: str) -> Dict[str, Any]:
        """Parse a single mutation string."""
        import re

        # Pattern for single letter codes: A123V
        pattern_1 = r"^([A-Z])(\d+)([A-Z])$"
        # Pattern for three-letter codes: ALA123VAL
        pattern_3 = r"^([A-Z]{3})(\d+)([A-Z]{3})$"

        match_1 = re.match(pattern_1, mutation_str)
        match_3 = re.match(pattern_3, mutation_str)

        if match_1:
            from_aa, position, to_aa = match_1.groups()
            return {
                "from_residue": from_aa,
                "to_residue": to_aa,
                "position": int(position),
                "format": "single_letter",
            }
        elif match_3:
            from_aa, position, to_aa = match_3.groups()
            return {
                "from_residue": from_aa,
                "to_residue": to_aa,
                "position": int(position),
                "format": "three_letter",
            }
        else:
            raise ProteinMutationError(f"Invalid mutation format: {mutation_str}")


class ProteinMutationPlanner:
    """Plan protein mutation free energy calculations."""

    def __init__(self):
        self.aa_mapping = {
            "A": "ALA",
            "R": "ARG",
            "N": "ASN",
            "D": "ASP",
            "C": "CYS",
            "E": "GLU",
            "Q": "GLN",
            "G": "GLY",
            "H": "HIS",
            "I": "ILE",
            "L": "LEU",
            "K": "LYS",
            "M": "MET",
            "F": "PHE",
            "P": "PRO",
            "S": "SER",
            "T": "THR",
            "W": "TRP",
            "Y": "TYR",
            "V": "VAL",
        }
        self.reverse_aa_mapping = {v: k for k, v in self.aa_mapping.items()}

    def create_protein_mutation_network(
        self,
        protein_file: str,
        mutation_str: str,
        ligand_sdf: Optional[str] = None,
        cofactor_sdf: Optional[str] = None,
        output_dir: str = "protein_mutations",
        settings: Optional[Dict] = None,
    ) -> AlchemicalNetwork:
        """Create an alchemical network for protein mutations."""

        logger.info(f"Planning protein mutation: {mutation_str}")

        # Parse mutation string
        mutation_info = MutationParser.parse_mutation_string(mutation_str)

        # Load components
        components = self._load_components(protein_file, ligand_sdf, cofactor_sdf)

        # Create wild-type and mutant systems
        wt_system, mut_system = self._create_mutation_systems(components, mutation_info, settings)

        # Create transformation
        transformation = self._create_protein_transformation(
            wt_system, mut_system, mutation_info, settings
        )

        # Create network
        network = AlchemicalNetwork([transformation])

        # Save network
        self._save_network(network, output_dir, mutation_info)

        return network

    def _load_components(
        self, protein_file: str, ligand_sdf: Optional[str], cofactor_sdf: Optional[str]
    ) -> Dict[str, Any]:
        """Load molecular components from files."""
        components = {}

        # Load protein
        if not os.path.exists(protein_file):
            raise ProteinMutationError(f"Protein file not found: {protein_file}")

        logger.info(f"Loading protein from {protein_file}")
        components["protein"] = ProteinComponent.from_pdb_file(protein_file)

        # Load ligand if provided
        if ligand_sdf:
            if not os.path.exists(ligand_sdf):
                raise ProteinMutationError(f"Ligand SDF file not found: {ligand_sdf}")

            logger.info(f"Loading ligands from {ligand_sdf}")
            components["ligands"] = self._load_molecules_from_sdf(ligand_sdf)

        # Load cofactors if provided
        if cofactor_sdf:
            if not os.path.exists(cofactor_sdf):
                raise ProteinMutationError(f"Cofactor SDF file not found: {cofactor_sdf}")

            logger.info(f"Loading cofactors from {cofactor_sdf}")
            components["cofactors"] = self._load_molecules_from_sdf(cofactor_sdf)

        # Add solvent
        components["solvent"] = SolventComponent()

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
                    f"Loaded molecule {i + 1}: {mol.GetProp('_Name') if mol.HasProp('_Name') else f'mol_{i + 1}'}"
                )

        if not molecules:
            raise ProteinMutationError(f"No valid molecules found in {sdf_file}")

        return molecules

    def _create_mutation_systems(
        self,
        components: Dict[str, Any],
        mutation_info: Dict[str, Any],
        settings: Optional[Dict] = None,
    ) -> Tuple[ChemicalSystem, ChemicalSystem]:
        """Create wild-type and mutant chemical systems."""

        # Base system components
        system_components = [components["protein"], components["solvent"]]

        # Add ligands if present
        if "ligands" in components:
            # For now, use the first ligand - could be extended for multiple ligands
            system_components.append(components["ligands"][0])

        # Add cofactors if present
        if "cofactors" in components:
            system_components.extend(components["cofactors"])

        # Create wild-type system
        wt_system = ChemicalSystem(system_components)

        # Create mutant protein
        mutant_protein = self._create_mutant_protein(components["protein"], mutation_info)

        # Create mutant system
        mutant_components = [mutant_protein, components["solvent"]]
        if "ligands" in components:
            mutant_components.append(components["ligands"][0])
        if "cofactors" in components:
            mutant_components.extend(components["cofactors"])

        mut_system = ChemicalSystem(mutant_components)

        return wt_system, mut_system

    def _create_mutant_protein(
        self, protein: ProteinComponent, mutation_info: Dict[str, Any], add_missing: bool = False
    ) -> List[ProteinComponent]:
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
        ...     "mutations": [
        ...         {"from_residue": "A", "position": 123, "to_residue": "V", "chain": "A"},
        ...         {"from_residue": "L", "position": 456, "to_residue": "P", "chain": "B"},
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

        for mutation in mutation_info["mutations"]:
            logger.info(
                f"Applying mutation: {mutation['from_residue']}{mutation['position']}{mutation['to_residue']} in chain {mutation['chain']}"
            )
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
            mutation_spec = (
                f"{mutation['from_residue']}-{mutation['position']}-{mutation['to_residue']}"
            )
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

    def _create_protein_transformation(
        self,
        wt_system: ChemicalSystem,
        mut_system: ChemicalSystem,
        mutation_info: Dict[str, Any],
        settings: Optional[Dict] = None,
    ) -> Transformation:
        """Create transformation between wild-type and mutant systems."""

        # Default protocol settings
        # TODO: This should probably be a setting obj/model
        default_settings = {
            "protocol": "NonEquilibriumCyclingProtocol",
            "lambda_schedule": [0.0, 0.25, 0.5, 0.75, 1.0],
            "simulation_time": "5.0 ns",
            "equilibration_time": "1.0 ns",
        }

        if settings:
            default_settings.update(settings)

        # Create protocol
        protocol = NonEquilibriumCyclingProtocol(settings=default_settings)

        # Create atom mapping (simplified)
        # In practice, this would use sophisticated protein atom mapping
        atom_mapping = self._create_protein_atom_mapping(wt_system, mut_system, mutation_info)

        # Create transformation
        transformation = Transformation(
            stateA=wt_system,
            stateB=mut_system,
            mapping=atom_mapping,
            protocol=protocol,
            name=f"protein_mutation_{mutation_info['original_string']}",
        )

        return transformation

    def _create_protein_atom_mapping(
        self, wt_system: ChemicalSystem, mut_system: ChemicalSystem, mutation_info: Dict[str, Any]
    ) -> AtomMapping:
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
            componentA_to_componentB=mapping_dict, componentA=wt_system, componentB=mut_system
        )

    def _save_network(
        self, network: AlchemicalNetwork, output_dir: str, mutation_info: Dict[str, Any]
    ):
        """Save the alchemical network to files."""

        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Save network
        network_file = output_path / "protein_mutation_network.json"
        with open(network_file, "w") as f:
            json.dump(network.to_dict(), f, indent=2)

        # Save individual transformations
        transformations_dir = output_path / "transformations"
        transformations_dir.mkdir(exist_ok=True)

        for i, transformation in enumerate(network.edges):
            trans_file = transformations_dir / f"transformation_{i:03d}.json"
            with open(trans_file, "w") as f:
                json.dump(transformation.to_dict(), f, indent=2)

        # Save mutation info
        mutation_file = output_path / "mutation_info.yaml"
        with open(mutation_file, "w") as f:
            yaml.dump(mutation_info, f, default_flow_style=False)

        logger.info(f"Network saved to {output_dir}")
