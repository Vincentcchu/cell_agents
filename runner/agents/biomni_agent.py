"""Biomni command wrapper."""
from __future__ import annotations

from ..registry import register_agent
from .command_agent import CommandAgent


@register_agent("biomni")
class BiomniAgent(CommandAgent):
    """Runs the Biomni agent using the configured command."""
