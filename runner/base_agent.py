"""Base abstractions used by the centralized agent runner."""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional


@dataclass(slots=True)
class AgentRequest:
    """Normalized payload passed to each agent wrapper."""

    data_path: Path
    task: str
    params: Dict[str, Any] = field(default_factory=dict)
    working_dir: Optional[Path] = None
    env: Dict[str, str] = field(default_factory=dict)

    def resolve_working_dir(self) -> Path:
        """Return a directory that exists on disk for the agent to use."""
        if self.working_dir is not None:
            self.working_dir.mkdir(parents=True, exist_ok=True)
            return self.working_dir
        default_dir = self.data_path.parent / "tmp"
        default_dir.mkdir(parents=True, exist_ok=True)
        return default_dir


@dataclass(slots=True)
class AgentResult:
    """Return object describing what an agent produced."""

    success: bool
    message: str = ""
    artifacts: List[Path] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)


class BaseAgent(ABC):
    """Shared interface that every concrete agent wrapper must implement."""

    name: str

    def __init__(self, name: str, config: Dict[str, Any] | None = None) -> None:
        self.name = name
        self.config = config or {}

    @abstractmethod
    def run(self, request: AgentRequest) -> AgentResult:  # pragma: no cover - interface only
        """Execute the agent using the normalized request payload."""

    def _merge_env(self, request: AgentRequest, extra_env: Optional[Dict[str, str]] = None) -> Dict[str, str]:
        """Combine runner-level and agent-level env overrides."""
        merged = dict(request.env)
        if extra_env:
            merged.update({k: v for k, v in extra_env.items() if v is not None})
        return merged
