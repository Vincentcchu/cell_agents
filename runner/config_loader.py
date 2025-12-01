"""Utility helpers for loading runner configuration files."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict

import yaml

RUNNER_ROOT = Path(__file__).resolve().parent
CONFIG_DIR = RUNNER_ROOT / "config"
DEFAULT_AGENT_CONFIG = CONFIG_DIR / "agent_configs.yaml"
DEFAULT_SECRET_FILE = CONFIG_DIR / "secrets.json"


def load_yaml(path: Path) -> Dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Missing configuration file: {path}")
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def load_json(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def load_agent_configs(path: Path | None = None) -> Dict[str, Any]:
    return load_yaml(path or DEFAULT_AGENT_CONFIG)


def load_secrets(path: Path | None = None) -> Dict[str, str]:
    data = load_json(path or DEFAULT_SECRET_FILE)
    if not isinstance(data, dict):
        raise ValueError("Secrets file must contain a JSON object")
    return {str(k): str(v) for k, v in data.items()}
