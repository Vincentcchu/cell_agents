"""Wrapper around the CASSIA Python API for batch inference."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CASSIA_PATH = PROJECT_ROOT / "agents" / "CASSIA" / "CASSIA_python"
if CASSIA_PATH.exists():
    sys.path.insert(0, str(CASSIA_PATH))

try:
    from CASSIA import runCASSIA_batch
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "Could not import CASSIA. Run `pip install -e agents/CASSIA/CASSIA_python` to install dependencies."
    ) from exc


def main() -> None:
    parser = argparse.ArgumentParser(description="Invoke CASSIA batch annotation")
    parser.add_argument("--data", required=True, help="Path to the marker CSV produced by preprocessing")
    parser.add_argument("--task", required=True, help="High-level description (for logging)")
    parser.add_argument("--output", default=None, help="Where to write the JSON output")
    parser.add_argument("--model", default=None)
    parser.add_argument("--provider", default=None)
    parser.add_argument("--tissue", default=None)
    parser.add_argument("--species", default=None)
    parser.add_argument("--n-genes", type=int, default=50)
    parser.add_argument("--temperature", type=float, default=0.0)
    parser.add_argument("--max-workers", type=int, default=6)
    parser.add_argument("--max-retries", type=int, default=1)
    parser.add_argument("--validator", default="v1", help="Validator involvement flag")
    args = parser.parse_args()

    marker_path = Path(args.data).resolve()
    if not marker_path.exists():
        raise SystemExit(f"Marker file not found: {marker_path}")

    output_file = args.output or "cassia_results.json"
    model = args.model or "google/gemini-2.5-flash-preview"
    provider = args.provider or "openrouter"
    tissue = args.tissue or "lung"
    species = args.species or "human"

    results = runCASSIA_batch(
        marker=str(marker_path),
        output_name=output_file,
        n_genes=args.n_genes,
        model=model,
        temperature=args.temperature,
        tissue=tissue,
        species=species,
        provider=provider,
        max_workers=args.max_workers,
        max_retries=args.max_retries,
        validator_involvement=args.validator,
    )

    output_path = Path(output_file)
    if not output_path.is_absolute():
        output_path = marker_path.parent / output_path
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(results, handle, indent=2)
    print(f"CASSIA results saved to {output_path}")


if __name__ == "__main__":
    main()
