"""Simple registry for agent wrappers."""
from __future__ import annotations

from typing import Dict, Iterable, Type

from .base_agent import BaseAgent

_AGENT_TYPES: Dict[str, Type[BaseAgent]] = {}


def register_agent(name: str):
    """Decorator used by agent wrapper classes."""

    def decorator(cls: Type[BaseAgent]) -> Type[BaseAgent]:
        key = name.lower()
        if key in _AGENT_TYPES:
            raise ValueError(f"Agent '{name}' already registered")
        _AGENT_TYPES[key] = cls
        return cls

    return decorator


def get_agent_class(name: str) -> Type[BaseAgent]:
    key = name.lower()
    if key not in _AGENT_TYPES:
        raise KeyError(f"Unknown agent '{name}'. Registered: {sorted(_AGENT_TYPES)}")
    return _AGENT_TYPES[key]


def list_agents() -> Iterable[str]:
    return sorted(_AGENT_TYPES)
