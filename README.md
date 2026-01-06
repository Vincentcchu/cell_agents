# LLM Cell-Type Annotation Benchmarking

This repository benchmarks multiple Large Language Model (LLM) agents on cell type annotation tasks using single-cell RNA-seq data. Each agent is independently developed and has its own preprocessing requirements, execution methods, and output formats.

## Repository Structure

```
cell_agents/
├── agents/              # Agent implementations (one folder per agent)
├── data/                # Shared datasets
├── preprocessing/       # Data preprocessing pipelines
├── outputs/             # Agent predictions and results
├── analysis/            # Evaluation notebooks and visualizations
├── runner/              # Optional unified runner interface
└── logs/                # Execution logs
```

## Quick Start

### 1. Data Location

All shared datasets are in the [`data/`](data/) folder:
- `dataset_restricted.h5ad` - Main restricted dataset
- `dataset_restricted_with_labels.h5ad` - Dataset with annotations
- See [data/README.md](data/README.md) for details

### 2. Running an Agent

Each agent has its own folder under [`agents/`](agents/). See the agent-specific README for instructions:

| Agent | Entrypoint | Documentation |
|-------|-----------|---------------|
| **BioMaster** | [`agents/biomaster/run.py`](agents/biomaster/run.py) | [README](agents/biomaster/README.md) / [Quick Start](agents/biomaster/QUICKSTART.md) |
| **Biomni** | [`agents/Biomni/run.ipynb`](agents/Biomni/run.ipynb) | [README](agents/Biomni/README.md) / [Quick Start](agents/Biomni/QUICKSTART.md) |
| **CASSIA** | [`agents/CASSIA/CASSIA_python/`](agents/CASSIA/CASSIA_python/) | [README](agents/CASSIA/README.md) / [Quick Start](agents/CASSIA/QUICKSTART.md) |
| **CellTypeAgent** | [`agents/CellTypeAgent/CellTypeAgent/`](agents/CellTypeAgent/CellTypeAgent/) | [README](agents/CellTypeAgent/README.md) / [Quick Start](agents/CellTypeAgent/QUICKSTART.md) |
| **mLLMCellType** | [`agents/mLLMCellType/run.py`](agents/mLLMCellType/run.py) | [Quick Start](agents/mLLMCellType/QUICKSTART.md) |

### 3. Preprocessing Data

Some agents require specific input formats. Preprocessing scripts are in [`preprocessing/`](preprocessing/):
- **Agent-specific**: [`preprocessing/agent_specific/`](preprocessing/agent_specific/)
- **Shared utilities**: [`preprocessing/shared/`](preprocessing/shared/)

See [preprocessing/README.md](preprocessing/README.md) for details.

### 4. Outputs

Agent predictions are saved to [`outputs/{agent_name}/`](outputs/):
- `outputs/biomaster/` - BioMaster results
- `outputs/biomni/` - Biomni results
- `outputs/cassia/` - CASSIA results
- `outputs/celltypeagent/` - CellTypeAgent results
- `outputs/mllmcelltype/` - mLLMCellType results

### 5. Evaluation

Evaluation and analysis tools are in [`analysis/`](analysis/):
- [`analysis/evaluator2.0.ipynb`](analysis/evaluator2.0.ipynb) - Main evaluation notebook
- [`analysis/umap_embeddings/`](analysis/umap_embeddings/) - UMAP visualizations
- See [analysis/README.md](analysis/README.md) for metrics and usage

## Workflow Overview

```
1. Prepare data → preprocessing/ (if needed for specific agent)
2. Run agent    → agents/{agent}/run.py or run.ipynb
3. Check output → outputs/{agent}/
4. Evaluate     → analysis/evaluator2.0.ipynb
```

## Design Principles

**Agent Independence**: Each agent is developed independently with its own:
- Preprocessing requirements
- Execution methods
- Dependencies
- Output formats

**No Forced Standardization**: Agents are not forced into a common interface. The optional `runner/` provides a convenience layer but is not required.

**Minimal Path Updates**: After reorganization, agents mostly reference data using relative paths from their own directories.

## Getting Started with a Specific Agent

1. **Choose an agent** from the table above
2. **Read its README** in `agents/{agent}/README.md`
3. **Prepare data** if preprocessing is required (see agent README)
4. **Run the agent** following instructions in its README
5. **Find outputs** in `outputs/{agent}/`
6. **Evaluate** using notebooks in `analysis/`

## Additional Resources

- **Agent overview**: [agents/README.md](agents/README.md)
- **Data documentation**: [data/README.md](data/README.md)
- **Preprocessing guide**: [preprocessing/README.md](preprocessing/README.md)
- **Analysis guide**: [analysis/README.md](analysis/README.md)
- **Unified runner** (optional): [runner/README.md](runner/README.md)

## Notes

- Agents may be git submodules pointing to external repositories
- Some agents have multiple experimental run notebooks (archived in `experiments/`)
- Logs are collected in `logs/` when using the unified runner
