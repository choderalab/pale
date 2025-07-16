"""
Current place to host custom ProtocolUnits that can be used for experimenting with running
specific units and parts of a gufe.ProtocolDAG.
"""

from typing import Any

from gufe import Context
from feflow.protocols.nonequilibrium_cycling import SetupUnit


class CustomSetupUnit(SetupUnit):
    """
    SetupUnit for feflow.protocol.NonEquilibriumCyclingProtocol that mimicks a SetupUnit for that
    protocol based on the contents of its shared context directory.
    """

    def _execute(ctx: Context, **inputs) -> dict[str, Any]:
        pass
