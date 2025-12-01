"""Simplified CLI that runs the mLLMCellType consensus annotation pipeline."""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import List

PROJECT_ROOT = Path(__file__).resolve().parents[2]
MLLM_ROOT = PROJECT_ROOT / "agents" / "mLLMCellType"
if MLLM_ROOT.exists():
    sys.path.insert(0, str(MLLM_ROOT))


def ensure_api_keys(json_path: Path | None) -> None:
    if json_path and json_path.exists():
        with json_path.open("r", encoding="utf-8") as handle:
            data = json.load(handle)
        for env_key, json_key in {
            "OPENAI_API_KEY": "openai_api_key",
            "OPENROUTER_API_KEY": "openrouter_api_key",
            "ANTHROPIC_API_KEY": "anthropic_api_key",
        }.items():
            if env_key not in os.environ and json_key in data:
                os.environ[env_key] = data[json_key]


def parse_models(raw: str) -> List[dict]:
    models: List[dict] = []
    if not raw:
        return models
    for chunk in raw.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        if ":" in chunk:
            provider, model = chunk.split(":", 1)
        else:
            provider, model = "openrouter", chunk
        models.append({"provider": provider, "model": model})
    return models


def compute_markers(adata, cluster_key: str, top_genes: int) -> dict:
    import scanpy as sc

    if cluster_key not in adata.obs:
        raise SystemExit(f"Cluster key '{cluster_key}' missing from AnnData object")
    sc.tl.rank_genes_groups(adata, cluster_key, method="wilcoxon")
    marker_genes = {}
    categories = list(adata.obs[cluster_key].cat.categories)
    for cluster in categories:
        df = sc.get.rank_genes_groups_df(adata, group=cluster).head(top_genes)
        marker_genes[str(cluster)] = df["names"].tolist()
    return marker_genes


def ensure_basic_processing(adata, cluster_key: str) -> None:
    import scanpy as sc

    if 'log1p' not in adata.uns:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    if 'X_pca' not in adata.obsm:
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        sc.pp.pca(adata, use_highly_variable=True)
    if 'neighbors' not in adata.uns:
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
    if cluster_key not in adata.obs:
        sc.tl.leiden(adata, resolution=0.8, key_added=cluster_key)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run mLLMCellType consensus annotation")
    parser.add_argument("--data", required=True, help="Path to the h5ad file")
    parser.add_argument("--task", required=True, help="Description logged alongside the run")
    parser.add_argument("--output", default=None)
    parser.add_argument("--cluster-key", default=None)
    parser.add_argument("--top-genes", type=int, default=10)
    parser.add_argument("--species", default=None)
    parser.add_argument("--tissue", default=None)
    parser.add_argument("--models", default=None)
    parser.add_argument("--consensus", type=float, default=0.7)
    parser.add_argument("--rounds", type=int, default=2)
    parser.add_argument("--api-config", default=None)
    args = parser.parse_args()

    try:
        import scanpy as sc
    except ImportError as exc:  # pragma: no cover
        raise SystemExit("scanpy is required for the mLLM runner") from exc

    try:
        from mllmcelltype import interactive_consensus_annotation
    except ImportError as exc:  # pragma: no cover
        raise SystemExit("Install the mllmcelltype package to use this runner") from exc

    data_path = Path(args.data).resolve()
    if not data_path.exists():
        raise SystemExit(f"h5ad file not found: {data_path}")

    output_path_value = args.output or "mllm_annotations.h5ad"
    cluster_key = args.cluster_key or "leiden"
    species = args.species or "human"
    tissue = args.tissue or "blood"
    api_config_path = Path(args.api_config) if args.api_config else (MLLM_ROOT / "api_keys.json")

    ensure_api_keys(api_config_path)

    adata = sc.read_h5ad(data_path)
    ensure_basic_processing(adata, cluster_key)
    marker_genes = compute_markers(adata, cluster_key, args.top_genes)

    models = parse_models(args.models or "")
    if not models:
        models = [
            {"provider": "openrouter", "model": "deepseek/deepseek-r1-0528:free"},
            {"provider": "openrouter", "model": "google/gemini-2.0-flash-exp:free"},
        ]

    results = interactive_consensus_annotation(
        marker_genes=marker_genes,
        species=species,
        tissue=tissue,
        models=models,
        consensus_threshold=args.consensus,
        max_discussion_rounds=args.rounds,
    )

    annotations = results["consensus"]
    adata.obs["mllm_cell_type"] = adata.obs[cluster_key].astype(str).map(annotations)
    adata.obs["mllm_confidence"] = adata.obs[cluster_key].astype(str).map(results.get("consensus_proportion", {}))

    output_path = Path(output_path_value)
    if not output_path.is_absolute():
        output_path = data_path.parent / output_path
    adata.write(output_path)
    print(f"mLLMCellType annotations written to {output_path}")


if __name__ == "__main__":
    main()
