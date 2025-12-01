"""Bridge script for running CellTypeAgent predictions on custom datasets."""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CELLTYPE_ROOT = PROJECT_ROOT / "agents" / "CellTypeAgent"
PACKAGE_ROOT = CELLTYPE_ROOT  # Outer folder contains the Python package
if PACKAGE_ROOT.exists():
    sys.path.insert(0, str(PACKAGE_ROOT))

try:
    from CellTypeAgent.get_prediction import get_prediction
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "Could not import CellTypeAgent. Run `pip install -e agents/CellTypeAgent` to expose the package."
    ) from exc


DATASETS_DIR = CELLTYPE_ROOT / "CellTypeAgent" / "data" / "GPTCellType" / "datasets"
DATASETS_DIR.mkdir(parents=True, exist_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run CellTypeAgent on a formatted CSV file")
    parser.add_argument("--data", required=True, help="Path to the formatted CSV generated during preprocessing")
    parser.add_argument("--task", required=True, help="Description logged alongside the run")
    parser.add_argument("--dataset-name", default=None, help="Dataset name to register inside CellTypeAgent")
    parser.add_argument("--model", default=None)
    parser.add_argument("--top-n", type=int, default=3, help="Number of candidates per cluster")
    parser.add_argument("--max-markers", type=int, default=None)
    parser.add_argument("--species", default=None)
    parser.add_argument("--log-dir", default=None)
    args = parser.parse_args()

    source = Path(args.data).resolve()
    if not source.exists():
        raise SystemExit(f"Input file not found: {source}")

    dataset_name = args.dataset_name or source.stem
    target = DATASETS_DIR / f"{dataset_name}.csv"
    shutil.copyfile(source, target)

    log_dir_value = args.log_dir or "runner_logs"
    log_dir = f"logs/{log_dir_value}"
    model = args.model or "o3-mini-2025-01-31"
    species = args.species or "human"
    data, output_csv = get_prediction(
        model=model,
        dataset_name=dataset_name,
        number_of_candidates=args.top_n,
        data_dir="data/GPTCellType/datasets",
        log_dir=log_dir,
        max_markers=args.max_markers,
        species=species,
    )

    print(f"CellTypeAgent predictions saved to {output_csv}")


if __name__ == "__main__":
    main()
