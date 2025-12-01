"""Invoke the BioMaster multi-agent pipeline with standardized arguments."""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
BIOMASTER_PATH = PROJECT_ROOT / "agents" / "biomaster"
if BIOMASTER_PATH.exists():
    sys.path.insert(0, str(BIOMASTER_PATH))


def env_or_default(name: str, default: str | None = None) -> str:
    value = os.environ.get(name)
    if value is None:
        if default is None:
            raise SystemExit(f"Missing required environment variable: {name}")
        return default
    return value


def str_to_bool(value: str | None, fallback: bool) -> bool:
    if value is None:
        return fallback
    return value.lower() in {"1", "true", "yes", "on"}


def main() -> None:
    parser = argparse.ArgumentParser(description="Invoke BioMaster via a unified CLI")
    parser.add_argument("--task", required=True, help="Goal prompt for BioMaster")
    parser.add_argument("--data", required=True, help="Path to the prepared dataset")
    parser.add_argument("--description", default="Preprocessed dataset ready for BioMaster")
    parser.add_argument("--skip-plan", action="store_true", help="Skip the PLAN step and execute tasks only")
    args = parser.parse_args()

    try:
        from agents.Biomaster import Biomaster
    except ImportError as exc:  # pragma: no cover
        raise SystemExit("Could not import BioMaster. Did you run `pip install -e agents/biomaster`?") from exc

    api_key = env_or_default("BIOMASTER_API_KEY")
    base_url = env_or_default("BIOMASTER_BASE_URL", "https://api.openai.com/v1")
    embedding_api_key = os.environ.get("BIOMASTER_EMBEDDING_API_KEY", api_key)
    embedding_base_url = os.environ.get("BIOMASTER_EMBEDDING_BASE_URL", base_url)
    model = os.environ.get("BIOMASTER_MODEL", "o3-mini-2025-01-31")
    tool_model = os.environ.get("BIOMASTER_TOOL_MODEL", model)
    embedding_model = os.environ.get("BIOMASTER_EMBEDDING_MODEL", "text-embedding-3-large")
    execute_debug = str_to_bool(os.environ.get("BIOMASTER_EXECUTOR"), False)
    use_ollama = str_to_bool(os.environ.get("BIOMASTER_USE_OLLAMA"), False)
    ollama_base_url = os.environ.get("BIOMASTER_OLLAMA_BASE_URL", "http://localhost:11434")
    run_id = os.environ.get("BIOMASTER_RUN_ID", "runner")
    output_dir = os.environ.get("BIOMASTER_OUTPUT_DIR", str(BIOMASTER_PATH / "output"))

    data_path = Path(args.data).resolve()
    datalist = [f"{data_path}: {args.description}"]

    manager = Biomaster(
        api_key=api_key,
        base_url=base_url,
        Model=model,
        tool_model=tool_model,
        embedding_model=embedding_model,
        embedding_base_url=embedding_base_url,
        embedding_api_key=embedding_api_key,
        excutor=execute_debug,
        id=run_id,
        output_dir=output_dir,
        use_ollama=use_ollama,
        ollama_base_url=ollama_base_url,
    )

    if not args.skip_plan:
        manager.execute_PLAN(args.task, datalist)
    results = manager.execute_TASK(datalist)
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
