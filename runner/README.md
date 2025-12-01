# Centralized Agent Runner

This package provides a thin abstraction layer that lets you trigger any of the agent sub-repositories from a single CLI entry point or via Python imports. Preprocessing pipelines remain separate; each agent wrapper expects data that has already been formatted for that agent.

## Planned layout

```
runner/
├── README.md                 # You are here
├── __init__.py               # Convenience exports
├── base_agent.py             # Shared interface + helpers
├── config/
│   ├── agent_configs.yaml    # Agent-specific settings
│   └── secrets.example.json  # Document expected API keys
├── agents/
│   ├── __init__.py
│   ├── biomni_agent.py
│   ├── biomaster_agent.py
│   ├── cassia_agent.py
│   ├── celltype_agent.py
│   └── mllm_agent.py
└── run_agent.py              # CLI for `python -m runner.run_agent`
```

Each wrapper implements a small CLI script under `runner/scripts/` and is invoked through the shared `CommandAgent` helper. Agent-specific keywords (e.g., tissue, model) are passed via the `--param key=value` flags on the runner CLI and are exposed to templates as `{param_key}` placeholders.

## Configuration & secrets

1. Copy `runner/config/secrets.example.json` to `runner/config/secrets.json` and fill in your API keys:
	 ```bash
	 cp runner/config/secrets.example.json runner/config/secrets.json
	 $EDITOR runner/config/secrets.json
	 ```
2. Review `runner/config/agent_configs.yaml` to tweak commands, defaults, and environment variables per agent. The `command` field accepts either a full string (parsed with `shlex.split`) or a YAML list of arguments.

Secrets referenced under `env.from_secrets` are injected as environment variables before each command is executed. Static entries (for example default models) can be placed under `env.static`.

## Usage

Run any agent through the shared CLI:

```bash
python -m runner.run_agent --list  # show available agents

python -m runner.run_agent \
	--agent biomni \
	--data data/dataset_restricted.h5ad \
	--task "Classify malignant vs non-malignant cells" \
	--param llm=claude-sonnet-4-20250514
```

Any `--param key=value` pair becomes available to the command template as `{param_key}`. Optional placeholders resolve to an empty string when not provided, so the underlying script should fall back to its internal defaults.
