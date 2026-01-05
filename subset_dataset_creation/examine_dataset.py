#!/usr/bin/env python3
"""
Script to examine the contents of integrated_with_quiescence.h5ad dataset
"""

import pandas as pd
import numpy as np
try:
    import anndata as ad
    import scanpy as sc
    ANNDATA_AVAILABLE = True
except ImportError:
    ANNDATA_AVAILABLE = False
    print("AnnData/Scanpy not available. Installing...")

def examine_h5ad_file(filepath):
    """Examine the contents of an h5ad file"""
    
    if not ANNDATA_AVAILABLE:
        print("Please install anndata and scanpy first:")
        print("pip install anndata scanpy")
        return
    
    try:
        # Load the data
        print(f"Loading {filepath}...")
        adata = ad.read_h5ad(filepath)
        
        print("\n" + "="*60)
        print("DATASET OVERVIEW")
        print("="*60)
        
        # Basic information
        print(f"Shape: {adata.shape[0]} cells × {adata.shape[1]} genes")
        print(f"Data type: {type(adata.X)}")
        if hasattr(adata.X, 'dtype'):
            print(f"Data dtype: {adata.X.dtype}")
        
        # Check if data is sparse
        from scipy import sparse
        if sparse.issparse(adata.X):
            print("Data format: Sparse matrix")
            print(f"Sparsity: {1 - adata.X.nnz / (adata.shape[0] * adata.shape[1]):.4f}")
        else:
            print("Data format: Dense matrix")
        
        print("\n" + "="*60)
        print("OBSERVATIONS (CELLS)")
        print("="*60)
        print(f"Number of cells: {adata.n_obs}")
        print(f"Cell metadata columns: {list(adata.obs.columns)}")
        
        # Show first few rows of cell metadata
        if not adata.obs.empty:
            print("\nFirst 5 rows of cell metadata:")
            print(adata.obs.head())
            
            # Show unique values for categorical columns
            print("\nUnique values in categorical columns:")
            for col in adata.obs.columns:
                if adata.obs[col].dtype == 'object' or adata.obs[col].dtype.name == 'category':
                    unique_vals = adata.obs[col].unique()
                    if len(unique_vals) <= 20:  # Only show if not too many unique values
                        print(f"  {col}: {unique_vals}")
                    else:
                        print(f"  {col}: {len(unique_vals)} unique values")
        
        print("\n" + "="*60)
        print("VARIABLES (GENES)")
        print("="*60)
        print(f"Number of genes: {adata.n_vars}")
        print(f"Gene metadata columns: {list(adata.var.columns)}")
        
        # Show first few rows of gene metadata
        if not adata.var.empty:
            print("\nFirst 5 rows of gene metadata:")
            print(adata.var.head())
        
        # Show gene names
        print(f"\nFirst 10 gene names: {adata.var_names[:10].tolist()}")
        
        print("\n" + "="*60)
        print("LAYERS")
        print("="*60)
        if adata.layers:
            print(f"Available layers: {list(adata.layers.keys())}")
            for layer_name in adata.layers.keys():
                layer_data = adata.layers[layer_name]
                if sparse.issparse(layer_data):
                    print(f"  {layer_name}: {layer_data.shape} (sparse)")
                else:
                    print(f"  {layer_name}: {layer_data.shape} (dense)")
        else:
            print("No additional layers found")
        
        print("\n" + "="*60)
        print("EMBEDDINGS & DIMENSIONALITY REDUCTIONS")
        print("="*60)
        if adata.obsm:
            print(f"Available embeddings: {list(adata.obsm.keys())}")
            for emb_name in adata.obsm.keys():
                emb_data = adata.obsm[emb_name]
                print(f"  {emb_name}: {emb_data.shape}")
        else:
            print("No embeddings found")
        
        print("\n" + "="*60)
        print("PAIRWISE ANNOTATIONS")
        print("="*60)
        if adata.obsp:
            print(f"Available pairwise annotations: {list(adata.obsp.keys())}")
            for obsp_name in adata.obsp.keys():
                obsp_data = adata.obsp[obsp_name]
                print(f"  {obsp_name}: {obsp_data.shape}")
        else:
            print("No pairwise annotations found")
        
        print("\n" + "="*60)
        print("UNSTRUCTURED ANNOTATIONS")
        print("="*60)
        if adata.uns:
            print(f"Available unstructured annotations: {list(adata.uns.keys())}")
            for uns_name in adata.uns.keys():
                uns_data = adata.uns[uns_name]
                print(f"  {uns_name}: {type(uns_data)}")
        else:
            print("No unstructured annotations found")
        
        # Data summary statistics
        print("\n" + "="*60)
        print("DATA STATISTICS")
        print("="*60)
        
        # Convert to dense for statistics if sparse
        if sparse.issparse(adata.X):
            # Sample a subset for statistics if dataset is large
            if adata.shape[0] > 1000:
                sample_idx = np.random.choice(adata.shape[0], 1000, replace=False)
                sample_data = adata.X[sample_idx, :].toarray()
                print("(Statistics computed on random sample of 1000 cells)")
            else:
                sample_data = adata.X.toarray()
        else:
            sample_data = adata.X
        
        print(f"Mean expression: {np.mean(sample_data):.4f}")
        print(f"Median expression: {np.median(sample_data):.4f}")
        print(f"Max expression: {np.max(sample_data):.4f}")
        print(f"Min expression: {np.min(sample_data):.4f}")
        print(f"Fraction of zeros: {np.mean(sample_data == 0):.4f}")
        
        # Check for specific quiescence-related information
        print("\n" + "="*60)
        print("QUIESCENCE-RELATED INFORMATION")
        print("="*60)
        
        # Look for quiescence-related columns
        quiescence_cols = []
        for col in adata.obs.columns:
            if 'quiesc' in col.lower() or 'dormant' in col.lower() or 'cycle' in col.lower() or 'phase' in col.lower():
                quiescence_cols.append(col)
        
        if quiescence_cols:
            print("Found quiescence-related columns:")
            for col in quiescence_cols:
                print(f"  {col}: {adata.obs[col].value_counts().to_dict()}")
        else:
            print("No obvious quiescence-related columns found in cell metadata")
            print("Checking for cell cycle or state-related columns...")
            
            state_cols = []
            for col in adata.obs.columns:
                if any(keyword in col.lower() for keyword in ['state', 'cluster', 'type', 'group', 'label']):
                    state_cols.append(col)
            
            if state_cols:
                print("Found potential cell state columns:")
                for col in state_cols[:5]:  # Show first 5 to avoid too much output
                    unique_vals = adata.obs[col].unique()
                    if len(unique_vals) <= 10:
                        print(f"  {col}: {unique_vals}")
                    else:
                        print(f"  {col}: {len(unique_vals)} unique values")
        
        return adata
        
    except Exception as e:
        print(f"Error reading file: {e}")
        return None

if __name__ == "__main__":
    filepath = "integrated_with_quiescence.h5ad"
    
    # Redirect output to file
    import sys
    with open('dataset_analysis.txt', 'w') as f:
        sys.stdout = f
        adata = examine_h5ad_file(filepath)
        sys.stdout = sys.__stdout__  # Reset stdout
    
    print("Analysis complete! Results saved to dataset_analysis.txt")
    print("You can now view the full analysis with: cat dataset_analysis.txt")
