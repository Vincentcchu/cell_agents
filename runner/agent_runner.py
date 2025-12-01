"""High-level entry point for executing agents via a unified interface."""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

from .base_agent import AgentRequest, AgentResult
from .config_loader import (
    CONFIG_DIR,
    DEFAULT_AGENT_CONFIG,
    DEFAULT_SECRET_FILE,
    load_agent_configs,
    load_secrets,
)
from .registry import get_agent_class


class AgentRunner:
    def __init__(
        self,
        config_path: Optional[Path] = None,
        secrets_path: Optional[Path] = None,
    ) -> None:
        self.repo_root = Path(__file__).resolve().parents[1]
        self.config_path = config_path or DEFAULT_AGENT_CONFIG
        self.secrets_path = secrets_path or DEFAULT_SECRET_FILE
        self.agent_configs = load_agent_configs(self.config_path)
        self.secrets = load_secrets(self.secrets_path)
        self._instances: Dict[str, Any] = {}

    def available_agents(self) -> Dict[str, Dict[str, Any]]:
        return self.agent_configs.get("agents", {})

    def get_agent(self, agent_name: str):
        key = agent_name.lower()
        if key in self._instances:
            return self._instances[key]
        agent_cfg = self._get_agent_config(key)
        agent_cls = get_agent_class(key)
        instance = agent_cls(name=key, config=agent_cfg)
        self._instances[key] = instance
        return instance

    def run(
        self,
        agent_name: str,
        data_path: str | Path,
        task: str,
        params: Optional[Dict[str, Any]] = None,
        working_dir: Optional[str | Path] = None,
    ) -> AgentResult:
        agent = self.get_agent(agent_name)
        cfg = self._get_agent_config(agent.name)
        payload = AgentRequest(
            data_path=self._resolve_path(data_path),
            task=task,
            params=params or {},
            working_dir=self._resolve_working_dir(cfg, working_dir),
            env=self._build_env(cfg),
        )
        return agent.run(payload)

    # Internal helpers -------------------------------------------------

    def _resolve_path(self, path_like: str | Path) -> Path:
        path = Path(path_like)
        if not path.is_absolute():
            path = (self.repo_root / path).resolve()
        if not path.exists():
            raise FileNotFoundError(f"Input data not found: {path}")
        return path

    def _resolve_working_dir(self, cfg: Dict[str, Any], override: Optional[str | Path]) -> Optional[Path]:
        if override is not None:
            return Path(override)
        if "working_dir" in cfg:
            return self._ensure_dir(cfg["working_dir"])
        return None

    def _ensure_dir(self, path_like: str) -> Path:
        path = Path(path_like)
        if not path.is_absolute():
            path = (CONFIG_DIR.parent / path).resolve()
        path.mkdir(parents=True, exist_ok=True)
        return path

    def _build_env(self, cfg: Dict[str, Any]) -> Dict[str, str]:
        env_cfg = cfg.get("env", {})
        env: Dict[str, str] = {}
        for key, value in env_cfg.get("static", {}).items():
            env[key] = str(value)
        for key, secret_name in env_cfg.get("from_secrets", {}).items():
            if secret_name not in self.secrets:
                raise KeyError(
                    f"Missing secret '{secret_name}' required for env var '{key}'."
                    f" Update {self.secrets_path} or adjust the agent config."
                )
            env[key] = self.secrets[secret_name]
        return env

    def _get_agent_config(self, agent_name: str) -> Dict[str, Any]:
        agents_section = self.agent_configs.get("agents", {})
        if agent_name not in agents_section:
            raise KeyError(
                f"Agent '{agent_name}' missing from {self.config_path}."
                f" Available entries: {sorted(agents_section)}"
            )
        return agents_section[agent_name]
