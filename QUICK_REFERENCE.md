# Quick Reference Guide

## Running Each Agent

### BioMaster
```bash
cd agents/biomaster
cp config.yaml.template config.yaml
# Edit config.yaml with your API keys
python run.py config.yaml
```

### Biomni
```bash
cd agents/Biomni
jupyter notebook run.ipynb
```

### CASSIA
```bash
# First, preprocess data:
cd preprocessing/agent_specific/cassia
jupyter notebook convert_format.ipynb

# Then run CASSIA:
cd ../../../agents/CASSIA/CASSIA_python/CASSIA
python run_full_test.py
```

### CellTypeAgent
```bash
# First, preprocess data:
cd preprocessing/agent_specific/celltypeagent
jupyter notebook format_convert_celltypeagent.ipynb

# Then run CellTypeAgent:
cd ../../../agents/CellTypeAgent/CellTypeAgent
# Follow CellTypeAgent-specific instructions
```

### mLLMCellType
```bash
cd agents/mLLMCellType
# Edit api_keys.json with your API keys
python run.py
```

## Directory Quick Reference

| What | Where |
|------|-------|
| Agent code | `agents/{agent}/` |
| Data files | `data/` |
| Agent-specific preprocessing | `preprocessing/agent_specific/{agent}/` |
| Shared preprocessing | `preprocessing/shared/{utility}/` |
| Agent outputs | `outputs/{agent}/` |
| Evaluation notebook | `analysis/evaluator2.0.ipynb` |
| UMAP visualizations | `analysis/umap_embeddings/` |
| Logs | `logs/` |

## Common Tasks

### Run evaluation:
```bash
cd analysis
jupyter notebook evaluator2.0.ipynb
```

### Remove labels from dataset:
```bash
cd preprocessing/shared/remove_label
jupyter notebook "remove label.ipynb"
```

### Create dataset subset:
```bash
cd preprocessing/shared/subset_dataset_creation
python create_debug_dataset.py --fraction 0.1
```

### Convert gene data to h5ad:
```bash
cd preprocessing/shared/gene_data_to_h5ad
jupyter notebook karaayvas2018_h5ad_curate.ipynb
```

## Path Conventions

From agent directories (`agents/{agent}/`):
- Data: `../../data/dataset_restricted.h5ad`
- Outputs: `../../outputs/{agent}/`

From preprocessing directories:
- Data: `../../data/`  (from `preprocessing/agent_specific/`)
- Data: `../../../data/` (from nested directories)

From analysis directory (`analysis/`):
- Data: `../data/`
- Outputs: `../outputs/{agent}/`

## Verification

Run the verification script to check the repository structure:
```bash
./verify_reorganization.sh
```

## Documentation

- Top-level overview: [README.md](README.md)
- Agent guides: `agents/{agent}/QUICKSTART.md`
- Data info: [data/README.md](data/README.md)
- Preprocessing info: [preprocessing/README.md](preprocessing/README.md)
- Analysis info: [analysis/README.md](analysis/README.md)
- Full reorganization details: [REORGANIZATION_SUMMARY.md](REORGANIZATION_SUMMARY.md)
