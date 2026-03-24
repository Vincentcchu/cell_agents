#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Iterable


VALID_SUBFOLDERS = ("clustered", "non_clustered")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Aggregate CELLxGENE CSV files within dataset subfolders while skipping metadata rows, "
            "and optionally filter by Tissue column."
        )
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("."),
        help="Root directory containing dataset folders (default: current directory).",
    )
    parser.add_argument(
        "--datasets",
        nargs="+",
        default=None,
        help=(
            "One or more dataset folder names to process (e.g., breast brain). "
            "If omitted, all dataset folders under --root are processed."
        ),
    )
    parser.add_argument(
        "--subfolders",
        nargs="+",
        choices=VALID_SUBFOLDERS,
        default=list(VALID_SUBFOLDERS),
        help="Which subfolders to process (default: clustered non_clustered).",
    )
    parser.add_argument(
        "--skip-rows",
        type=int,
        default=9,
        help="Number of metadata rows to skip at top of each CSV (default: 9).",
    )
    parser.add_argument(
        "--tissue-filter",
        nargs="+",
        default=None,
        help=(
            "Filter values for the 'Tissue' column applied after aggregation. "
            "Example: --tissue-filter 'breast' 'adipose tissue'"
        ),
    )
    parser.add_argument(
        "--case-sensitive",
        action="store_true",
        help="Make --tissue-filter matching case-sensitive (default: case-insensitive).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("aggregated_output"),
        help="Directory where aggregated CSV files are written (default: aggregated_output).",
    )
    return parser.parse_args()


def discover_datasets(root: Path, requested: list[str] | None) -> list[Path]:
    if requested:
        datasets = [root / name for name in requested]
    else:
        datasets = [path for path in root.iterdir() if path.is_dir()]

    valid = [path for path in datasets if any((path / sf).is_dir() for sf in VALID_SUBFOLDERS)]
    return sorted(valid)


def get_csv_files(folder: Path) -> list[Path]:
    return sorted(path for path in folder.glob("*.csv") if path.is_file())


def read_header(csv_path: Path, skip_rows: int) -> list[str] | None:
    with csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
        for _ in range(skip_rows):
            next(handle, None)
        reader = csv.reader(handle)
        return next(reader, None)


def build_unified_header(csv_files: Iterable[Path], skip_rows: int) -> list[str]:
    unified: list[str] = []
    seen = set()

    for csv_path in csv_files:
        header = read_header(csv_path, skip_rows)
        if not header:
            continue
        for column in header:
            if column not in seen:
                seen.add(column)
                unified.append(column)

    return unified


def iterate_rows(csv_path: Path, skip_rows: int) -> Iterable[dict[str, str]]:
    with csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
        for _ in range(skip_rows):
            next(handle, None)
        reader = csv.DictReader(handle)
        for row in reader:
            if row is None:
                continue
            yield row


def aggregate_files(csv_files: list[Path], output_csv: Path, skip_rows: int) -> tuple[int, list[str]]:
    header = build_unified_header(csv_files, skip_rows)
    if not header:
        return 0, []

    total_rows = 0
    with output_csv.open("w", encoding="utf-8", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, fieldnames=header)
        writer.writeheader()

        for csv_path in csv_files:
            for row in iterate_rows(csv_path, skip_rows):
                normalized = {column: row.get(column, "") for column in header}
                writer.writerow(normalized)
                total_rows += 1

    return total_rows, header


def filter_by_tissue(
    source_csv: Path,
    target_csv: Path,
    tissues: list[str],
    case_sensitive: bool,
) -> int:
    if case_sensitive:
        allowed = set(tissues)

        def is_allowed(value: str) -> bool:
            return value in allowed

    else:
        allowed = {value.lower() for value in tissues}

        def is_allowed(value: str) -> bool:
            return value.lower() in allowed

    matched_rows = 0
    with source_csv.open("r", encoding="utf-8", newline="") as in_handle, target_csv.open(
        "w", encoding="utf-8", newline=""
    ) as out_handle:
        reader = csv.DictReader(in_handle)
        if not reader.fieldnames:
            return 0

        writer = csv.DictWriter(out_handle, fieldnames=reader.fieldnames)
        writer.writeheader()

        for row in reader:
            tissue_value = (row.get("Tissue") or "").strip()
            if is_allowed(tissue_value):
                writer.writerow(row)
                matched_rows += 1

    return matched_rows


def process_subfolder(
    dataset_dir: Path,
    subfolder: str,
    skip_rows: int,
    output_dir: Path,
    tissue_filter: list[str] | None,
    case_sensitive: bool,
) -> None:
    source_dir = dataset_dir / subfolder
    csv_files = get_csv_files(source_dir)
    if not csv_files:
        print(f"[SKIP] {source_dir}: no CSV files found")
        return

    output_dir.mkdir(parents=True, exist_ok=True)

    base_name = f"{dataset_dir.name}_{subfolder}_aggregated"
    aggregate_path = output_dir / f"{base_name}.csv"
    row_count, header = aggregate_files(csv_files, aggregate_path, skip_rows)

    if not header:
        print(f"[SKIP] {source_dir}: no tabular header found after skipping {skip_rows} rows")
        return

    print(f"[OK] {source_dir}: merged {len(csv_files)} file(s), {row_count} rows -> {aggregate_path}")

    if tissue_filter:
        filtered_path = output_dir / f"{base_name}_filtered.csv"
        matched = filter_by_tissue(
            source_csv=aggregate_path,
            target_csv=filtered_path,
            tissues=tissue_filter,
            case_sensitive=case_sensitive,
        )
        print(
            f"[OK] {source_dir}: filtered by Tissue ({', '.join(tissue_filter)}) -> "
            f"{matched} rows in {filtered_path}"
        )


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    output_dir = args.output_dir.resolve()

    datasets = discover_datasets(root, args.datasets)
    if not datasets:
        print(f"No dataset folders found under {root}")
        return

    for dataset_dir in datasets:
        for subfolder in args.subfolders:
            if (dataset_dir / subfolder).is_dir():
                process_subfolder(
                    dataset_dir=dataset_dir,
                    subfolder=subfolder,
                    skip_rows=args.skip_rows,
                    output_dir=output_dir,
                    tissue_filter=args.tissue_filter,
                    case_sensitive=args.case_sensitive,
                )


if __name__ == "__main__":
    main()


# python aggregate_cellxgene_csvs.py --datasets brain breast --subfolders non_clustered --output-dir aggregated_output
# python aggregate_cellxgene_csvs.py --datasets breast --tissue-filter "abdomen"