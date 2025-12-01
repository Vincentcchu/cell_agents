"""Lightweight wrapper to invoke Biomni's Python API."""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
BIOMNI_PATH = PROJECT_ROOT / "agents" / "Biomni"
if BIOMNI_PATH.exists():
    sys.path.insert(0, str(BIOMNI_PATH))


def build_prompt(task: str, data_path: Path) -> str:
    return f"{task}\nDataset location: {data_path}"


def main() -> None:
    parser = argparse.ArgumentParser(description="Invoke Biomni using the unified runner")
    parser.add_argument("--task", required=True, help="Goal to hand off to Biomni")
    parser.add_argument("--data", required=True, help="Path to the prepared data")
    parser.add_argument("--llm", default=None)
    args = parser.parse_args()

    try:
        from biomni.agent import A1
    except ImportError as exc:  # pragma: no cover - import side effect
        raise SystemExit(
            "Could not import biomni.agent. Ensure the Biomni submodule is installed (pip install -e agents/Biomni)."
        ) from exc

    llm = args.llm or os.environ.get("BIOMNI_LLM", "claude-sonnet-4-20250514")
    data_path = Path(args.data).resolve()
    data_dir = data_path if data_path.is_dir() else data_path.parent
    prompt = build_prompt(args.task, data_path)

    agent = A1(path=str(data_dir), llm=llm)
    agent.go(prompt)


if __name__ == "__main__":
    main()
