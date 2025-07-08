#!/usr/bin/env python3
"""
Script to verify that cell type annotations were properly removed from the dataset
"""

import sys
import pandas as pd

try:
    import anndata as ad
    ANNDATA_AVAILABLE = True
except ImportError:
    ANNDATA_AVAILABLE = False
    print("AnnData not available. Please install first.")
    sys.exit(1)

def verify_annotations_removed(original_filepath, modified_filepath):
    """
    Compare the original and modified datasets to verify cell type annotations were removed
    """
    print(f"Loading original dataset: {original_filepath}")
    original_adata = ad.read_h5ad(original_filepath)
    
    print(f"Loading modified dataset: {modified_filepath}")
    modified_adata = ad.read_h5ad(modified_filepath)
    
    print("\n--- COMPARISON SUMMARY ---")
    print(f"Original dataset shape: {original_adata.shape[0]} cells × {original_adata.shape[1]} genes")
    print(f"Modified dataset shape: {modified_adata.shape[0]} cells × {modified_adata.shape[1]} genes")
    
    print("\n--- METADATA COLUMNS ---")
    original_cols = set(original_adata.obs.columns)
    modified_cols = set(modified_adata.obs.columns)
    
    removed_cols = original_cols - modified_cols
    
    print(f"Original columns: {len(original_cols)}")
    print(f"Modified columns: {len(modified_cols)}")
    print(f"Removed columns: {len(removed_cols)}")
    
    if removed_cols:
        print("\nColumns that were removed:")
        for col in sorted(removed_cols):
            print(f"  - {col}")
    else:
        print("\nNo columns were removed!")
    
    # Check for information in .uns
    if 'removed_annotations' in modified_adata.uns:
        print("\n--- REMOVAL METADATA ---")
        removed_info = modified_adata.uns['removed_annotations']
        print(f"Columns removed according to metadata: {removed_info.get('columns', [])}")
        print(f"Unstructured annotations removed: {removed_info.get('uns_keys', [])}")
        print(f"Date of removal: {removed_info.get('date', 'unknown')}")
    
    return {
        'original_shape': original_adata.shape,
        'modified_shape': modified_adata.shape,
        'removed_columns': removed_cols
    }

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Verify cell type annotations were removed")
    parser.add_argument("--original", "-o", required=True,
                       help="Original h5ad file path")
    parser.add_argument("--modified", "-m", required=True,
                       help="Modified h5ad file path")
    
    args = parser.parse_args()
    
    verify_annotations_removed(args.original, args.modified)
