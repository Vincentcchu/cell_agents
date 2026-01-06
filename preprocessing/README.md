# Data Preprocessing

This directory contains preprocessing pipelines for preparing data in formats required by different agents.

## Structure

```
preprocessing/
├── README.md                    # This file
├── agent_specific/              # Agent-specific format conversions
│   ├── cassia/                  # CASSIA preprocessing
│   │   └── convert_format.ipynb
│   └── celltypeagent/           # CellTypeAgent preprocessing
│       └── format_convert_celltypeagent.ipynb
└── shared/                      # Shared preprocessing utilities
    ├── remove_label/            # Remove cell type labels
    ├── subset_dataset_creation/ # Create dataset subsets
    └── gene_data_to_h5ad/       # Convert gene data to h5ad
```

## Agent-Specific Preprocessing

Some agents require data in specific formats. These notebooks handle the conversion:

### CASSIA Format Conversion

**Location**: [`agent_specific/cassia/convert_format.ipynb`](agent_specific/cassia/convert_format.ipynb)

Converts standard h5ad files to CASSIA-specific format.

**Usage**:
```bash
cd agent_specific/cassia/
jupyter notebook convert_format.ipynb
```

**Input**: `../../data/dataset_restricted.h5ad`  
**Output**: CASSIA-formatted files for agent input

### CellTypeAgent Format Conversion

**Location**: [`agent_specific/celltypeagent/format_convert_celltypeagent.ipynb`](agent_specific/celltypeagent/format_convert_celltypeagent.ipynb)

Converts standard h5ad files to CellTypeAgent-specific format.

**Usage**:
```bash
cd agent_specific/celltypeagent/
jupyter notebook format_convert_celltypeagent.ipynb
```

**Input**: `../../data/dataset_restricted.h5ad`  
**Output**: CellTypeAgent-formatted files

## Shared Preprocessing Utilities

These utilities work with any agent and provide common data manipulation tasks.

### 1. Remove Cell Type Labels

**Location**: [`shared/remove_label/`](shared/remove_label/)

Removes cell type annotations from annotated datasets to create blind test sets.

**Documentation**: [shared/remove_label/README.md](shared/remove_label/README.md)

**Usage**:
```bash
cd shared/remove_label/
jupyter notebook "remove label.ipynb"
```

**Use case**: Create test datasets where ground truth is hidden from agents

### 2. Create Dataset Subsets

**Location**: [`shared/subset_dataset_creation/`](shared/subset_dataset_creation/)

Creates smaller subsets of large datasets for testing and debugging.

**Documentation**: [shared/subset_dataset_creation/README_toy_dataset.md](shared/subset_dataset_creation/README_toy_dataset.md)

**Usage**:
```bash
cd shared/subset_dataset_creation/
python create_debug_dataset.py
```

**Features**:
- Random sampling
- Stratified sampling (preserves cell type distributions)
- Configurable sample size
- Quality checks

### 3. Gene Data to h5ad

**Location**: [`shared/gene_data_to_h5ad/`](shared/gene_data_to_h5ad/)

Converts raw gene expression data (MTX format) to AnnData h5ad format.

**Usage**:
```bash
cd shared/gene_data_to_h5ad/
jupyter notebook karaayvas2018_h5ad_curate.ipynb
```

**Input**: 
- `Exp_data_TPM.mtx` - Expression matrix
- `Genes.txt` - Gene names

**Output**: Standard h5ad file

## Preprocessing Workflow

### For new datasets:

1. **Convert raw data to h5ad** (if needed):
   ```bash
   cd shared/gene_data_to_h5ad/
   # Edit notebook to point to your data
   jupyter notebook karaayvas2018_h5ad_curate.ipynb
   ```

2. **Remove labels** (if creating test set):
   ```bash
   cd shared/remove_label/
   # Configure input/output paths in notebook
   jupyter notebook "remove label.ipynb"
   ```

3. **Create subset** (optional, for testing):
   ```bash
   cd shared/subset_dataset_creation/
   python create_debug_dataset.py --fraction 0.1
   ```

4. **Agent-specific conversion** (if needed):
   ```bash
   # For CASSIA:
   cd agent_specific/cassia/
   jupyter notebook convert_format.ipynb
   
   # For CellTypeAgent:
   cd agent_specific/celltypeagent/
   jupyter notebook format_convert_celltypeagent.ipynb
   ```

## Which Agents Need Preprocessing?

| Agent | Preprocessing Required | Script Location |
|-------|----------------------|-----------------|
| BioMaster | ❌ No | Uses standard h5ad |
| Biomni | ❌ No | Uses standard h5ad |
| CASSIA | ✅ Yes | `agent_specific/cassia/` |
| CellTypeAgent | ✅ Yes | `agent_specific/celltypeagent/` |
| mLLMCellType | ❌ No | Uses standard h5ad |

## Adding New Preprocessing

To add preprocessing for a new agent:

1. Create directory: `agent_specific/{agent_name}/`
2. Add conversion script/notebook
3. Document inputs, outputs, and usage
4. Update this README

## Notes

- All preprocessing scripts assume data is in `../../data/`
- Converted data should be used directly by agents, not saved back to `data/`
- Preprocessing is idempotent - running multiple times produces same results
- Always verify data integrity after conversion
