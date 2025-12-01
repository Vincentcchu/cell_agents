# Test Data for Machine Learning Models

This directory contains datasets prepared for testing machine learning models on single-cell RNA-seq data.

## Files

- **dataset_debug.h5ad**: Original dataset with cell type annotations
- **dataset_debug_no_celltype.h5ad**: Modified dataset with cell type annotations removed
- **dataset_modification_summary.md**: Detailed summary of modifications made to remove cell type information

## Purpose

The modified dataset (`dataset_debug_no_celltype.h5ad`) is designed for ML model testing by removing all cell type annotations while preserving the underlying expression data. This creates a scenario where:

1. The ground truth labels are not available to the model
2. The expression matrix and technical metadata are preserved
3. You can evaluate model performance by comparing predictions to the original dataset labels

## Usage

### For Model Testing

```python
import scanpy as sc

# Load the unlabeled dataset
adata = sc.read_h5ad('dataset_debug_no_celltype.h5ad')

# Run your ML model
# ...

# Compare with ground truth (using original dataset separately)
ground_truth = sc.read_h5ad('dataset_debug.h5ad')
```

### Scripts Used

The following scripts were used to create and analyze these datasets:

- `remove_cell_type_annotations.py`: Removes cell type annotations from h5ad files
- `verify_annotations_removed.py`: Verifies that annotations were correctly removed
- `create_summary_document.py`: Creates a detailed summary of the modification process

## Ground Truth Access

For evaluation purposes, you can use the original `dataset_debug.h5ad` file, which contains the ground truth cell type annotations. This allows you to:

1. Train models on the unlabeled data
2. Generate predictions
3. Evaluate performance using the ground truth labels from the original dataset
