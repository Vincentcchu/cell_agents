"""Command-based agent wrappers."""
from __future__ import annotations

import os
import shlex
import subprocess
from pathlib import Path
from typing import Any, Dict, Iterable, List

from ..base_agent import AgentRequest, AgentResult, BaseAgent


class _SafeFormatDict(dict):
    def __missing__(self, key):  # pragma: no cover - simple fallback
        return ""

PROJECT_ROOT = Path(__file__).resolve().parents[2]


class CommandAgent(BaseAgent):
    """Agent wrapper that shells out to a configurable command."""

    def run(self, request: AgentRequest) -> AgentResult:
        command = self._build_command(request)
        cwd = self._resolve_cwd()
        working_dir = request.resolve_working_dir()
        env = os.environ.copy()
        env.update(request.env)
        env["AGENT_WORKDIR"] = str(working_dir)
        try:
            completed = subprocess.run(
                command,
                cwd=cwd,
                env=env,
                check=False,
            )
        except FileNotFoundError as exc:
            return AgentResult(success=False, message=f"Executable not found: {exc}")
        success = completed.returncode == 0
        message = f"Command exited with code {completed.returncode}"
        return AgentResult(
            success=success,
            message=message,
            metadata={
                "command": command,
                "cwd": str(cwd),
                "returncode": completed.returncode,
            },
        )

    # Helpers ---------------------------------------------------------

    def _build_command(self, request: AgentRequest) -> List[str]:
        template = self.config.get("command")
        if not template:
            raise ValueError(f"Agent '{self.name}' is missing a 'command' entry in the config")
        flat_context = self._build_context(request)
        safe_context = _SafeFormatDict(flat_context)
        if isinstance(template, str):
            rendered = template.format_map(safe_context)
            return shlex.split(rendered)
        return [str(token).format_map(safe_context) for token in template]

    def _build_context(self, request: AgentRequest) -> Dict[str, str]:
        context: Dict[str, str] = {
            "task": request.task,
            "data_path": str(request.data_path),
            "data_dir": str(request.data_path.parent),
            "working_dir": str(request.resolve_working_dir()),
            "agent_name": self.name,
        }
        params = request.params or {}
        for key, value in params.items():
            context[f"param_{key}"] = str(value)
        return context

    def _resolve_cwd(self) -> Path:
        cwd_config = self.config.get("cwd")
        if not cwd_config:
            return PROJECT_ROOT
        path = Path(cwd_config)
        if not path.is_absolute():
            path = (PROJECT_ROOT / path).resolve()
        return path
