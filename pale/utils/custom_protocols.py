"""
File to experiment with custom protocols that can be of utility for extending runs of previously
run protocols or to do some things different compared to "standard" OpenFE protocols.
"""

from os import PathLike
from typing import Optional

from feflow.protocols.nonequilibrium_cycling import NonEquilibriumCyclingProtocol, ResultUnit
from gufe import Settings, ChemicalSystem, ProtocolDAGResult, ProtocolUnit

from pale.utils import SetupFromShared


class NEqCyclingFromSetup(NonEquilibriumCyclingProtocol):
    """
    Class that inherits from the NEq Cycling protocol in order to create a new protocol that
    runs cycles from a existing successful setup unit result in a directory.

    Note that you need a fully successful run of a feflow.protocols.nonequilibrium_cycling.SetupUnit
    result in a directory for this protocol to work.
    """

    def __init__(self, settings: Settings):
        super().__init__(settings)

    def _create(
        self,
        stateA: ChemicalSystem,
        stateB: ChemicalSystem,
        shared_dir: Optional[PathLike] = None,
        extends: Optional[ProtocolDAGResult] = None,
    ) -> list[ProtocolUnit]:
        """
        Create a Protocol DAG from a custom setup from shared directory unit.
        """
        if shared_dir is None:
            raise ValueError("`shared_dir` is required for this Protocol. Please specify one.")

        setup = SetupFromShared(
            protocol=self,
            shared_dir=shared_dir,
            name="setup_from_shared",
        )

        num_cycles = self.settings.num_cycles

        # Create a cycling unit per cycle
        simulations = [
            self._simulation_unit(protocol=self, setup=setup, name=f"cycle_{replicate}")
            for replicate in range(num_cycles)
        ]

        end = ResultUnit(name="result", simulations=simulations)

        return [setup, *simulations, end]
