"""mLLMCellType command wrapper."""
from __future__ import annotations

from ..registry import register_agent
from .command_agent import CommandAgent


@register_agent("mllmcelltype")
class MLLMCellTypeAgent(CommandAgent):
    """Runs the mLLMCellType agent using the configured command."""
