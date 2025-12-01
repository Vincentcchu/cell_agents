"""CLI entry point for the centralized agent runner."""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict

from . import AgentRunner
from .agents import *  # noqa: F401,F403  # Ensure agent modules register themselves


def parse_params(pairs: list[str]) -> Dict[str, str]:
    parsed: Dict[str, str] = {}
    for item in pairs:
        if "=" not in item:
            raise ValueError(f"Invalid --param '{item}'. Expected KEY=VALUE format.")
        key, value = item.split("=", 1)
        parsed[key] = value
    return parsed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run an agent from the centralized registry")
    parser.add_argument("--agent", required=False, help="Agent name to invoke (e.g. biomni)")
    parser.add_argument("--data", help="Path to the preprocessed data file for the agent")
    parser.add_argument("--task", help="Natural-language task or goal for the agent")
    parser.add_argument(
        "--param",
        action="append",
        default=[],
        help="Additional key=value parameters forwarded to the agent",
    )
    parser.add_argument("--working-dir", help="Optional working directory override", default=None)
    parser.add_argument("--config", help="Custom agent config file", default=None)
    parser.add_argument("--secrets", help="Custom secrets JSON file", default=None)
    parser.add_argument(
        "--list",
        action="store_true",
        help="List available agents and exit",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    runner = AgentRunner(
        config_path=Path(args.config) if args.config else None,
        secrets_path=Path(args.secrets) if args.secrets else None,
    )

    if args.list or not args.agent:
        agents = runner.available_agents()
        lines = ["Registered agents:"]
        for name, cfg in agents.items():
            description = cfg.get("description", "")
            lines.append(f"- {name}: {description}")
        print("\n".join(lines))
        return

    params = parse_params(args.param)
    result = runner.run(
        agent_name=args.agent,
        data_path=args.data,
        task=args.task,
        params=params,
        working_dir=args.working_dir,
    )
    status = "succeeded" if result.success else "failed"
    print(f"Agent run {status}: {result.message}")


if __name__ == "__main__":
    main()
