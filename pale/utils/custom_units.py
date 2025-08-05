"""
Current place to host custom ProtocolUnits that can be used for experimenting with running
specific units and parts of a gufe.ProtocolDAG.
"""

from pathlib import Path
from typing import Any

from gufe import Context
from feflow.protocols.nonequilibrium_cycling import SetupUnit


class SetupFromShared(SetupUnit):
    """
    SetupUnit for feflow.protocol.NonEquilibriumCyclingProtocol that mimics a SetupUnit for that
    protocol based on the contents of its shared context directory.
    """

    def _execute(self, ctx: Context, shared_dir: Path, **inputs) -> dict[str, Any]:
        """
        The main idea here is to read a shared directory from a SetupUnit and restore the unit's
        outputs attribute accordingly such that it can be read by the CycleUnit.
        """
        import numpy as np

        system_outfile = shared_dir / "system.xml.bz2"
        state_outfile = shared_dir / "state.xml.bz2"
        integrator_outfile = shared_dir / "integrator.xml.bz2"
        htf_outfile = shared_dir / "hybrid_topology_factory.pickle"

        # We need to deserialize the HTF and extract the atom indices
        htf = np.load(htf_outfile)
        initial_atom_indices = htf.initial_atom_indices
        final_atom_indices = htf.final_atom_indices

        # It has to contain the same items as the base SetupUnit._execute
        outputs = {
            "system": system_outfile,
            "state": state_outfile,
            "integrator": integrator_outfile,
            "initial_atom_indices": initial_atom_indices,
            "final_atom_indices": final_atom_indices,
            "topology_path": htf_outfile,
        }

        return outputs
