"""CellTypeAgent command wrapper."""
from __future__ import annotations

from ..registry import register_agent
from .command_agent import CommandAgent


@register_agent("celltypeagent")
class CellTypeAgent(CommandAgent):
    """Runs CellTypeAgent using the configured command."""
