# Toy Dataset Creator for Single-Cell Data

This script creates a smaller version of the `integrated_with_quiescence.h5ad` dataset by randomly sampling a fraction of the cells. This is useful for testing code, developing new methods, or sharing a representative subset of the data.

## Features

- **Random sampling**: Choose a random subset of cells
- **Stratified sampling**: Preserves the distribution of important cell categories (e.g., disease status, cell types)
- **Configurable sample size**: Set any fraction of the original dataset
- **Quality checks**: Verifies the distribution of key categories
- **Metadata preservation**: All original metadata is preserved for sampled cells

## Usage

### Basic Usage (Stratified Sampling)

```bash
python create_toy_dataset.py
```

This will create a toy dataset with default settings:
- Takes 10% of the cells
- Stratifies based on 'disease' column
- Saves as `integrated_with_quiescence_toy_stratified.h5ad`
- Uses random seed 42 for reproducibility

### Random Sampling

```bash
python create_toy_dataset.py --random
```

This uses simple random sampling instead of stratified sampling.

### Custom Options

```bash
python create_toy_dataset.py --fraction 0.05 --stratify cell_type --output custom_toy_dataset.h5ad
```

This will:
- Take 5% of the cells
- Stratify based on 'cell_type' column
- Save to the specified output file

## All Options

```
--input, -i      Input h5ad file path
--output, -o     Output h5ad file path
--fraction, -f   Sampling fraction (default: 0.1)
--seed, -s       Random seed (default: 42)
--stratify, -t   Column name for stratified sampling (default: disease)
--random, -r     Use random sampling instead of stratified sampling
```

## How It Works

The script performs these steps:
1. Loads the original h5ad file
2. Either:
   - Randomly samples a fraction of cells, or
   - Performs stratified sampling to preserve category distributions
3. Compares the original and sampled distributions
4. Saves the result to a new h5ad file

## Example

```bash
# Create a toy dataset with 15% of cells, stratified by QuiescenceStatus
python create_toy_dataset.py -f 0.15 -t QuiescenceStatus -o quiescence_subset.h5ad
```

This will create a smaller dataset that preserves the proportion of quiescent, proliferating, and slow-cycling cells.
