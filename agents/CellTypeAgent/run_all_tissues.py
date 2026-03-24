#!/usr/bin/env python3
"""
Batch runner for CellTypeAgent over tissue-organized datasets.

Expected input structure:
    data/
      <tissue>/
        celltypeagent_format/
          <dataset>_formatted.csv

Outputs are organized as:
    <output_dir>/<tissue>/<dataset>/prediction/<model>/<timestamp>/...

Each run also writes a batch summary CSV/JSON in <output_dir>.
"""

import argparse
import csv
import json
import sys
import time
import traceback
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

CELLTYPE_AGENT_DIR = SCRIPT_DIR / "CellTypeAgent"
if str(CELLTYPE_AGENT_DIR) not in sys.path:
    sys.path.insert(0, str(CELLTYPE_AGENT_DIR))

from CellTypeAgent.get_prediction import get_prediction

DEFAULT_DATA_DIR_CANDIDATES = [
    SCRIPT_DIR / "data",
    SCRIPT_DIR.parents[2] / "data",
    SCRIPT_DIR.parents[3] / "data",
]
DEFAULT_DATA_DIR = next((path for path in DEFAULT_DATA_DIR_CANDIDATES if path.exists()), SCRIPT_DIR.parents[2] / "data")
DEFAULT_OUTPUT_DIR = SCRIPT_DIR / "batch_outputs"
DEFAULT_CSV_SUBDIR = "celltypeagent_format"


def discover_formatted_datasets(data_dir: Path, csv_subdir: str, file_pattern: str = "*.csv") -> List[Tuple[str, Path, str]]:
    datasets: List[Tuple[str, Path, str]] = []
    skip_dirs = {"__pycache__", ".git", ".ipynb_checkpoints"}

    for tissue_dir in sorted(data_dir.iterdir()):
        if not tissue_dir.is_dir() or tissue_dir.name in skip_dirs:
            continue

        candidate_dir = tissue_dir / csv_subdir
        if not candidate_dir.exists() or not candidate_dir.is_dir():
            continue

        for csv_file in sorted(candidate_dir.glob(file_pattern)):
            dataset_name = csv_file.stem
            if dataset_name.endswith("_all_datasets_formatted"):
                continue
            datasets.append((tissue_dir.name, csv_file, dataset_name))

    return datasets


def parse_start_from(value: Optional[str]) -> Tuple[Optional[str], Optional[str]]:
    if not value:
        return None, None
    if "/" in value:
        tissue, dataset = value.split("/", 1)
        return tissue.strip(), dataset.strip()
    return None, value.strip()


def apply_start_from(
    datasets: List[Tuple[str, Path, str]],
    start_from: Optional[str],
) -> List[Tuple[str, Path, str]]:
    if not start_from:
        return datasets

    start_tissue, start_dataset = parse_start_from(start_from)
    filtered: List[Tuple[str, Path, str]] = []
    found = False

    for tissue, csv_path, dataset_name in datasets:
        if not found:
            tissue_match = start_tissue is None or tissue == start_tissue
            dataset_match = (
                dataset_name == start_dataset
                or dataset_name.lower() == (start_dataset or "").lower()
                or (start_dataset or "").lower() in dataset_name.lower()
            )
            if tissue_match and dataset_match:
                found = True
                filtered.append((tissue, csv_path, dataset_name))
        else:
            filtered.append((tissue, csv_path, dataset_name))

    if not found:
        raise ValueError(f"Could not find start point: {start_from}")

    return filtered


def has_existing_run(output_dir: Path, tissue: str, dataset_name: str, model: str) -> bool:
    model_dir = output_dir / tissue / dataset_name / "prediction" / model
    if not model_dir.exists():
        return False
    return any(model_dir.glob("*/run_metrics.json"))


def write_batch_csv(results: List[Dict[str, Any]], csv_path: Path) -> None:
    if not results:
        return

    field_order = [
        "tissue",
        "dataset",
        "input_csv",
        "status",
        "success",
        "started_at",
        "finished_at",
        "wall_time_seconds",
        "output_prediction_csv",
        "run_metrics_json",
        "input_tokens",
        "output_tokens",
        "total_tokens",
        "preprocessing_time_seconds",
        "annotation_time_seconds",
        "total_time_seconds",
        "input_cost_usd",
        "output_cost_usd",
        "total_cost_usd",
        "error",
    ]

    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(csv_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=field_order)
        writer.writeheader()
        for row in results:
            writer.writerow({k: row.get(k, "") for k in field_order})


def main() -> None:
    parser = argparse.ArgumentParser(description="Run CellTypeAgent prediction across all tissue datasets")
    parser.add_argument("--data-dir", type=str, default=str(DEFAULT_DATA_DIR), help="Base data directory containing tissue subfolders")
    parser.add_argument("--csv-subdir", type=str, default=DEFAULT_CSV_SUBDIR, help="Subdirectory under each tissue containing CellTypeAgent formatted CSVs")
    parser.add_argument("--file-pattern", type=str, default="*.csv", help="Glob pattern used to pick dataset CSV files")
    parser.add_argument("--output-dir", type=str, default=str(DEFAULT_OUTPUT_DIR), help="Where outputs and batch summaries are written")

    parser.add_argument("--model", type=str, default="gpt-5.1", help="Model name passed to CellTypeAgent")
    parser.add_argument("--species", type=str, default="human", help="Species override for prompt generation")
    parser.add_argument("--top-n", type=int, default=3, help="Number of candidate cell types")
    parser.add_argument("--max-markers", type=int, default=None, help="Maximum markers per row (default: all)")
    parser.add_argument("--poor-performance-model", action="store_true", help="Enable strict count reminder in prompt")

    parser.add_argument("--input-cost-per-million", type=float, default=1.1, help="USD per 1M input tokens")
    parser.add_argument("--output-cost-per-million", type=float, default=4.4, help="USD per 1M output tokens")

    parser.add_argument("--tissues", nargs="+", default=None, help="Only process selected tissues")
    parser.add_argument("--start-from", type=str, default=None, help="Start from dataset (format: tissue/dataset or dataset)")
    parser.add_argument("--skip-existing", action="store_true", help="Skip datasets that already have run metrics for the selected model")
    parser.add_argument("--dry-run", action="store_true", help="Preview what would be run")

    args = parser.parse_args()

    data_dir = Path(args.data_dir).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if not data_dir.exists():
        raise FileNotFoundError(f"Data directory not found: {data_dir}")

    print(f"[INFO] Scanning: {data_dir} (subdir={args.csv_subdir}, pattern={args.file_pattern})")
    datasets = discover_formatted_datasets(data_dir, args.csv_subdir, args.file_pattern)

    if args.tissues:
        tissue_set = set(args.tissues)
        datasets = [entry for entry in datasets if entry[0] in tissue_set]
        print(f"[INFO] Tissue filter applied: {sorted(tissue_set)}")

    datasets = apply_start_from(datasets, args.start_from)

    if not datasets:
        print("[WARN] No datasets found with current filters")
        return

    print(f"[INFO] {len(datasets)} dataset(s) discovered")
    if args.dry_run:
        print("[INFO] Dry-run mode: no model calls will be made")

    results: List[Dict[str, Any]] = []
    batch_started = datetime.now().isoformat(timespec="seconds")

    for idx, (tissue, csv_path, dataset_name) in enumerate(datasets, start=1):
        print("\n" + "=" * 80)
        print(f"[{idx}/{len(datasets)}] {tissue}/{dataset_name}")
        print("=" * 80)

        started_at = datetime.now().isoformat(timespec="seconds")
        wall_start = time.perf_counter()

        result_row: Dict[str, Any] = {
            "tissue": tissue,
            "dataset": dataset_name,
            "input_csv": str(csv_path),
            "status": "pending",
            "success": False,
            "started_at": started_at,
            "finished_at": None,
            "wall_time_seconds": None,
            "output_prediction_csv": None,
            "run_metrics_json": None,
            "input_tokens": None,
            "output_tokens": None,
            "total_tokens": None,
            "preprocessing_time_seconds": None,
            "annotation_time_seconds": None,
            "total_time_seconds": None,
            "input_cost_usd": None,
            "output_cost_usd": None,
            "total_cost_usd": None,
            "error": None,
        }

        if args.skip_existing and has_existing_run(output_dir, tissue, dataset_name, args.model):
            result_row["status"] = "skipped_existing"
            result_row["success"] = True
            result_row["finished_at"] = datetime.now().isoformat(timespec="seconds")
            result_row["wall_time_seconds"] = round(time.perf_counter() - wall_start, 4)
            results.append(result_row)
            print("[SKIP] Existing run_metrics.json found for this dataset/model")
            continue

        if args.dry_run:
            result_row["status"] = "dry_run"
            result_row["success"] = True
            result_row["finished_at"] = datetime.now().isoformat(timespec="seconds")
            result_row["wall_time_seconds"] = round(time.perf_counter() - wall_start, 4)
            results.append(result_row)
            print(f"[DRY-RUN] Would run get_prediction for {csv_path}")
            continue

        try:
            data, prediction_csv_path = get_prediction(
                model=args.model,
                dataset_name=dataset_name,
                number_of_candidates=args.top_n,
                data_dir=str(csv_path.parent),
                log_dir=str(output_dir / tissue),
                poor_performance_model=args.poor_performance_model,
                max_markers=args.max_markers,
                species=args.species,
                input_cost_per_million=args.input_cost_per_million,
                output_cost_per_million=args.output_cost_per_million,
            )

            prediction_csv_path = Path(prediction_csv_path)
            run_dir = prediction_csv_path.parent
            metrics_path = run_dir / "run_metrics.json"

            dataset_root = output_dir / tissue / dataset_name
            dataset_root.mkdir(parents=True, exist_ok=True)
            top_txt = dataset_root / f"top_{args.top_n}_max_{args.max_markers}.txt"
            data["cell_type_pred"].to_csv(top_txt, index=False, header=False)

            result_row["status"] = "completed"
            result_row["success"] = True
            result_row["output_prediction_csv"] = str(prediction_csv_path)
            result_row["run_metrics_json"] = str(metrics_path) if metrics_path.exists() else None

            if metrics_path.exists():
                with open(metrics_path, "r", encoding="utf-8") as handle:
                    metrics = json.load(handle)

                result_row["input_tokens"] = metrics.get("input_tokens")
                result_row["output_tokens"] = metrics.get("output_tokens")
                result_row["total_tokens"] = metrics.get("total_tokens")
                result_row["preprocessing_time_seconds"] = metrics.get("preprocessing_time_seconds")
                result_row["annotation_time_seconds"] = metrics.get("annotation_time_seconds")
                result_row["total_time_seconds"] = metrics.get("total_time_seconds")
                result_row["input_cost_usd"] = metrics.get("input_cost_usd")
                result_row["output_cost_usd"] = metrics.get("output_cost_usd")
                result_row["total_cost_usd"] = metrics.get("total_cost_usd")

            print(f"[OK] Prediction CSV: {prediction_csv_path}")
            print(f"[OK] Dataset output folder: {dataset_root}")

        except Exception as exc:
            result_row["status"] = "failed"
            result_row["success"] = False
            result_row["error"] = f"{type(exc).__name__}: {exc}"
            print(f"[ERROR] {result_row['error']}")
            traceback.print_exc()

        result_row["finished_at"] = datetime.now().isoformat(timespec="seconds")
        result_row["wall_time_seconds"] = round(time.perf_counter() - wall_start, 4)
        results.append(result_row)

    batch_finished = datetime.now().isoformat(timespec="seconds")
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    summary = {
        "batch_start": batch_started,
        "batch_end": batch_finished,
        "model": args.model,
        "species": args.species,
        "input_cost_per_million": args.input_cost_per_million,
        "output_cost_per_million": args.output_cost_per_million,
        "total_runs": len(results),
        "successful": sum(1 for row in results if row.get("success")),
        "failed": sum(1 for row in results if row.get("status") == "failed"),
        "skipped": sum(1 for row in results if row.get("status") == "skipped_existing"),
        "dry_run": sum(1 for row in results if row.get("status") == "dry_run"),
        "results": results,
    }

    summary_json = output_dir / f"batch_summary_{timestamp}.json"
    summary_csv = output_dir / f"batch_summary_{timestamp}.csv"
    latest_json = output_dir / "BATCH_SUMMARY.json"
    latest_csv = output_dir / "BATCH_SUMMARY.csv"

    with open(summary_json, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    with open(latest_json, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    write_batch_csv(results, summary_csv)
    write_batch_csv(results, latest_csv)

    print("\n" + "=" * 80)
    print("BATCH COMPLETE")
    print("=" * 80)
    print(f"Total: {summary['total_runs']} | Successful: {summary['successful']} | Failed: {summary['failed']} | Skipped: {summary['skipped']}")
    print(f"Summary JSON: {summary_json}")
    print(f"Summary CSV : {summary_csv}")
    print(f"Latest JSON : {latest_json}")
    print(f"Latest CSV  : {latest_csv}")


if __name__ == "__main__":
    main()


# python run_all_tissues.py --model gpt-5.1 --species human
# Examples:
# python run_all_tissues.py --tissues brain breast
# python run_all_tissues.py --skip-existing
# python run_all_tissues.py --start-from brain/Data_Choudhury2022_Brain_formatted