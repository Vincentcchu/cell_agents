# Analysis & Evaluation

This directory contains evaluation notebooks, visualization tools, and analysis results for comparing agent performance.

## Contents

### Main Evaluation Notebook

**[`evaluator2.0.ipynb`](evaluator2.0.ipynb)**

Primary notebook for evaluating and comparing agent predictions against ground truth.

**Features**:
- Binary classification metrics (accuracy, precision, recall, F1)
- Multi-class classification evaluation
- Confusion matrix generation
- UMAP visualizations
- Cross-agent comparison

**Usage**:
```bash
jupyter notebook evaluator2.0.ipynb
```

**Inputs**:
- Agent predictions from `../outputs/{agent}/`
- Ground truth from `../data/dataset_restricted_with_labels.h5ad`

**Outputs**:
- Classification metrics
- Confusion matrices
- UMAP plots saved to `umap_embeddings/`
- Comparison CSV: `classification_for_evaluation.csv`

### UMAP Embeddings

**[`umap_embeddings/`](umap_embeddings/)**

Stores UMAP coordinates for visualization of agent predictions and ground truth.

**Files**:
- `{agent}_dataset_pred.npy` - Predictions from each agent
- `{agent}_dataset_truth.npy` - Corresponding ground truth labels
- `{agent}_dataset_cell_names.npy` - Cell identifiers

**Agents included**:
- `biomaster_*` - BioMaster results
- `biomni_*`, `biomni2_*` - Biomni results
- `mllmcelltype_*` - mLLMCellType results
- `dormancy_*`, `sarcoma_*` - Specialized analyses

### Classification Results

**[`classification_for_evaluation.csv`](classification_for_evaluation.csv)**

Summary CSV containing metrics for all agents, used for cross-comparison and reporting.

**Columns** (typical):
- Agent name
- Accuracy
- Precision
- Recall
- F1 score
- Additional task-specific metrics

### Token Counter

**[`token_counter.txt`](token_counter.txt)**

Tracks API token usage for LLM-based agents (cost estimation).

## Evaluation Workflow

### 1. Run agents and collect predictions:

```bash
# Predictions automatically saved to ../outputs/{agent}/
```

### 2. Run evaluation:

```bash
jupyter notebook evaluator2.0.ipynb
```

### 3. Review results:

- Check printed metrics in notebook
- View UMAP visualizations
- Compare across agents in `classification_for_evaluation.csv`

## Evaluation Metrics

### Binary Classification

For malignant vs. non-malignant classification:
- **Accuracy**: Overall correctness
- **Precision**: Positive predictive value
- **Recall**: Sensitivity/True positive rate
- **F1 Score**: Harmonic mean of precision and recall
- **Confusion Matrix**: True/False positives/negatives

### Multi-class Classification

For detailed cell type annotation:
- **Macro-averaged metrics**: Equal weight per class
- **Weighted metrics**: Weighted by class frequency
- **Per-class metrics**: Individual performance per cell type
- **Confusion matrix**: Full class-by-class breakdown

### Visualization

- **UMAP plots**: Compare predicted vs. true labels in low-dimensional space
- **Confusion matrices**: Heatmaps showing classification patterns
- **Distribution plots**: Class balance and prediction distributions

## UMAP Computation

UMAP embeddings are computed using:
```python
import umap

reducer = umap.UMAP(random_state=42, n_neighbors=15, min_dist=0.1)
embedding = reducer.fit_transform(expression_data)
```

**Configuration**:
- `random_state=42` - For reproducibility
- `n_neighbors=15` - Local structure preservation
- `min_dist=0.1` - Minimum distance between points

**Recomputation**:
Set `compute_umap=True` in evaluation functions to recompute. Otherwise, loads from `umap_embeddings/`.

## Adding New Analyses

To add evaluation for a new agent:

1. **Ensure output format**:
   - Agent saves predictions to `../outputs/{agent}/`
   - Predictions include cell IDs and predicted labels

2. **Load in evaluation notebook**:
   ```python
   import scanpy as sc
   pred_data = sc.read_h5ad('../outputs/{agent}/predictions.h5ad')
   ```

3. **Run evaluation functions**:
   ```python
   metrics = evaluate_binary_classification(
       pred_data, ground_truth,
       prediction_col='predicted_cell_type',
       truth_col='cell_type'
   )
   ```

4. **Save UMAP embeddings**:
   ```python
   np.save('umap_embeddings/{agent}_dataset_pred.npy', predictions)
   np.save('umap_embeddings/{agent}_dataset_truth.npy', ground_truth)
   ```

## Cross-Agent Comparison

The evaluation notebook includes functions for comparing multiple agents:

```python
# Compare multiple agents
agents = ['biomaster', 'biomni', 'cassia', 'celltypeagent', 'mllmcelltype']
compare_agents(agents, ground_truth, metric='f1_score')
```

Results are aggregated in `classification_for_evaluation.csv`.

## Reproducibility

To ensure reproducible evaluations:
- Use fixed random seeds (`random_state=42`)
- Save UMAP embeddings for consistent visualization
- Document preprocessing steps
- Version control evaluation notebooks

## Notes

- UMAP embeddings are large; consider `.gitignore` for binary files
- Token counting helps estimate API costs for LLM agents
- Evaluation notebook can be customized per analysis task
- Ground truth must be available for evaluation (not for agent input)
