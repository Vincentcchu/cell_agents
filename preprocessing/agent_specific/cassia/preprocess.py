import scanpy as sc
import CASSIA

INPUT_H5AD = "your_unlabelled_breast.h5ad"
OUTPUT_H5AD = "your_unlabelled_breast__prepped_for_cassia.h5ad"

# CASSIA marker extraction hyperparam (not in your snippet)
N_MARKERS_PER_CLUSTER = 50  # you can change; CASSIA examples commonly use ~50

def main():
    adata = sc.read_h5ad(INPUT_H5AD)

    # -------------------------
    # Basic preprocessing + clustering (only if missing)
    # -------------------------
    if "leiden" not in adata.obs.columns:
        # Normalize/log only if not already logged
        if "log1p" not in adata.uns:
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)

        # HVG + PCA only if missing PCA
        if "X_pca" not in adata.obsm:
            sc.pp.highly_variable_genes(
                adata,
                min_mean=0.0125,
                max_mean=3,
                min_disp=0.5
            )
            sc.pp.pca(adata, use_highly_variable=True)

        # Neighbors + Leiden
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
        sc.tl.leiden(adata, resolution=0.8, key_added="leiden")

    # -------------------------
    # Differential expression for marker genes (required by CASSIA bridge)
    # -------------------------
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")

    # -------------------------
    # Convert Scanpy markers -> CASSIA marker dataframe
    # (adds pct.1/pct.2 etc.; used by runCASSIA_* functions)
    # -------------------------
    markers_df = CASSIA.enhance_scanpy_markers(
        adata,
        cluster_col="leiden",
        n_genes=N_MARKERS_PER_CLUSTER
    )

    # Store markers inside the h5ad for easy notebook handoff
    adata.uns["CASSIA_markers"] = markers_df

    # Save
    adata.write_h5ad(OUTPUT_H5AD)
    print(f"Saved: {OUTPUT_H5AD}")
    print("Stored markers at: adata.uns['CASSIA_markers']")
    print("Cluster column: adata.obs['leiden']")

if __name__ == "__main__":
    main()
