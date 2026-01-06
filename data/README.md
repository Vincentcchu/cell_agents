# Data Directory

This directory contains shared datasets used for benchmarking cell type annotation agents.

## Available Datasets

### Primary Datasets

- **`dataset_restricted.h5ad`**
  - Main dataset with cell type labels removed (for agent prediction)
  - Format: AnnData h5ad
  - Used by: All agents (with preprocessing as needed)

- **`dataset_restricted_with_labels.h5ad`**
  - Same dataset with ground truth cell type annotations
  - Format: AnnData h5ad
  - Used for: Evaluation and comparison

### Processed Datasets

- **`integrated_with_quiescence.h5ad`**
  - Integrated dataset including quiescence markers
  
- **`integrated_with_quiescence_mapped.h5ad`**
  - Integrated dataset with cell type mappings

## Data Structure

Standard h5ad files contain:
- **`X`**: Gene expression matrix (cells × genes)
- **`obs`**: Cell-level metadata (cell IDs, annotations when present)
- **`var`**: Gene-level metadata (gene names, symbols)
- **`obsm`**: Multi-dimensional annotations (e.g., PCA, UMAP)
- **`uns`**: Unstructured metadata

## Usage

### Loading data in Python:

```python
import scanpy as sc

# Load dataset without labels (for prediction)
adata = sc.read_h5ad('data/dataset_restricted.h5ad')

# Load ground truth (for evaluation)
adata_truth = sc.read_h5ad('data/dataset_restricted_with_labels.h5ad')
```

### Data preprocessing:

Some agents require specific input formats. See:
- **Preprocessing scripts**: `../preprocessing/agent_specific/`
- **Shared utilities**: `../preprocessing/shared/`

## Data Preparation Pipeline

For creating test datasets or subsets:

1. **Remove labels**: See `../preprocessing/shared/remove_label/`
2. **Create subsets**: See `../preprocessing/shared/subset_dataset_creation/`
3. **Convert formats**: See `../preprocessing/shared/gene_data_to_h5ad/`

## Important Notes

- **Labels removed**: `dataset_restricted.h5ad` has cell type annotations removed for blind prediction
- **Ground truth**: Use `dataset_restricted_with_labels.h5ad` only for evaluation, not as agent input
- **Agent-specific formats**: Some agents need data conversion (see preprocessing documentation)

## Data Size

Check file sizes:
```bash
ls -lh *.h5ad
```

Typical sizes:
- Main datasets: 100MB - 1GB+
- Debug/toy datasets: 1-50MB

## Adding New Datasets

1. Place h5ad files in this directory
2. Update this README with dataset description
3. Update preprocessing scripts if format conversion is needed
4. Document in agent-specific README if applicable

## Privacy & Ethics

- Ensure datasets comply with data sharing agreements
- Remove or anonymize sensitive metadata
- Document data sources and citations appropriately
