#!/usr/bin/env python3
"""
Script to create a very small toy version of the h5ad dataset for testing purposes
"""

import pandas as pd
import numpy as np
import sys
import os
import scanpy as sc
import anndata as ad
from scipy import sparse

def create_tiny_toy_dataset(input_filepath, output_filepath=None, n_cells=1000, random_seed=42):
    """
    Create a tiny toy version with just a few cells for testing
    """
    try:
        np.random.seed(random_seed)
        
        # Default output path
        if output_filepath is None:
            output_filepath = "dataset_debug.h5ad"
            
        print(f"Loading dataset from {input_filepath}...")
        print("(This is only loading a small subset for quick testing)")
        
        # Use scanpy to read just enough cells
        # This avoids loading the entire dataset into memory
        adata = sc.read_h5ad(input_filepath, backed='r')
        
        # Get shape info without loading all data
        total_cells = adata.shape[0]
        n_cells = min(n_cells, total_cells)
        
        print(f"Creating tiny test dataset with {n_cells} cells (full dataset: {total_cells} cells)")
        
        # Sample indices
        indices = np.random.choice(total_cells, n_cells, replace=False)
        indices.sort()
        
        print("Using alternative approach to avoid backed mode issues...")
        
        # Close the backed file and reopen without backing
        adata.file.close()
        
        print("Loading a smaller portion of the data directly...")
        # Load the file directly but with as little as needed
        small_adata = sc.read_h5ad(input_filepath, backed=False, 
                               obs_names=indices)
        
        # Force to dense if very small and if it's a sparse matrix
        if n_cells <= 1000:
            try:
                if sparse.issparse(small_adata.X):
                    print("Converting to dense matrix format for small dataset...")
                    small_adata.X = small_adata.X.toarray()
            except Exception as e:
                print(f"Note: Could not convert to dense format: {e}")
                print("Continuing with original format.")
            
        # Add sampling information
        small_adata.uns['tiny_dataset_info'] = {
            'original_size': total_cells,
            'sampled_size': n_cells,
            'random_seed': random_seed,
            'purpose': 'testing'
        }
        
        # Save the tiny dataset
        print(f"Saving tiny test dataset to {output_filepath}...")
        small_adata.write_h5ad(output_filepath)
        
        print(f"Tiny test dataset created with {n_cells} cells")
        print(f"File saved to: {output_filepath}")
        
        return True
        
    except Exception as e:
        print(f"Error creating tiny dataset: {e}")
        return False

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Create a tiny test h5ad dataset")
    parser.add_argument("--input", "-i", default="integrated_with_quiescence.h5ad",
                       help="Input h5ad file path")
    parser.add_argument("--output", "-o", default="tiny_test_dataset.h5ad",
                       help="Output h5ad file path")
    parser.add_argument("--cells", "-c", type=int, default=1000,
                       help="Number of cells to include (default: 1000)")
    parser.add_argument("--seed", "-s", type=int, default=42,
                       help="Random seed for reproducibility")
    
    args = parser.parse_args()
    
    create_tiny_toy_dataset(
        input_filepath=args.input,
        output_filepath=args.output,
        n_cells=args.cells,
        random_seed=args.seed
    )
