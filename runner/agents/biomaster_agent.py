"""BioMaster command wrapper."""
from __future__ import annotations

from ..registry import register_agent
from .command_agent import CommandAgent


@register_agent("biomaster")
class BioMasterAgent(CommandAgent):
    """Runs the BioMaster agent via its configured command."""
