"""CASSIA command wrapper."""
from __future__ import annotations

from ..registry import register_agent
from .command_agent import CommandAgent


@register_agent("cassia")
class CassiaAgent(CommandAgent):
    """Runs the CASSIA agent using the configured command."""
