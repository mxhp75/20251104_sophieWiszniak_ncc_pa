#!/usr/bin/env python3
"""

=================

This script takes a filtered and normalised Seurat object, converted previously to a .h5ad file, and generates the paga graph (coarse-grained connectivity map between cell type clusters) and uses this information generate the FDG projections

Author: Melanie Smith
Email:  melanie.smith@adelaide.edu.au
Date:   2026-02-17

Usage:
    python this_script.py [options] [arguments]

Example:
    python this_script.py --input data.csv --output results/

Notes:
    - Main libraries used: ...
    - Python version: >= 3.9 recommended
"""

#-----------------------------------------------------------------------------------
#   Global configs
#-----------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------
#   Import core libraries
#-----------------------------------------------------------------------------------

import pandas as pd
import scanpy as sc # single-cell analysis toolkit
import matplotlib.pyplot as plt
import os

#-----------------------------------------------------------------------------------
#   Params
#-----------------------------------------------------------------------------------
fPCA=50 # set no. of pcs - this is used to populate the adata.obsm['X_pca'] field (run once)

# set the sample ID
# sample_id = "e11_control"
# sample_id = "e12_control"
# sample_id = "e11_ko"
sample_id = "e12_ko"

base_dir = "/home/melanie-smith/workDir/sophieWiszniak/20251104_sophieWiszniak_ncc_pa/outDir/09-rna_velocity"
# set the sample id and output directory
out_dir = "/home/melanie-smith/workDir/sophieWiszniak/20251104_sophieWiszniak_ncc_pa/outDir/09-rna_velocity"

# set the count level (ie gene or spliced/unspliced)
#count_level = "gene"
count_level = "spliced_unspliced"

#-----------------------------------------------------------------------------------
# create embeddings for dynamic velocity modelling
# PAGA informed Force Directed Graph
#-----------------------------------------------------------------------------------

# add the adata file - this is the processed and filtered RNA count matrix from Seurat converted to .h5ad
adata = sc.read_h5ad(os.path.join(base_dir, sample_id, f"{sample_id}_{count_level}.h5ad"))
print(adata)
print(adata.layers.keys())       # at least 'counts', possibly others
print(adata.obsm.keys())         # X_pca, X_umap, X_umap_spliced

# Compute PCA once with enough components (50) to cover all tested n_pcs values
sc.tl.pca(adata, svd_solver='arpack', n_comps=fPCA)

#-----------------------------------------------------------------------------------
# Explore the nPC and nNeighbours range
# PAGA informed Force Directed Graph
# Save output for all permutations
#-----------------------------------------------------------------------------------

pc_range = [10, 15, 20, 30] # this is how many PCs will be considered -> does not create a new adata.obsm['X_pca']
neighbour_range = [15, 25, 30, 50]

# Create output directory if it doesn't exist
save_dir = os.path.join(out_dir, sample_id, f"scanpy_output/celltypes")
os.makedirs(save_dir, exist_ok=True)

for n_pcs in pc_range:
    for n_neighbours in neighbour_range:
        adata_test = adata.copy()

        sc.pp.neighbors(adata_test, n_pcs=n_pcs, n_neighbors=n_neighbours)

        # Compute PAGA
        sc.tl.paga(adata_test, groups="new_celltypes")

        # Plot & SAVE PAGA
        sc.pl.paga(
            adata_test,
            threshold=0.03,
            color="new_celltypes",
            node_size_scale=4,
            fontsize=5,
            title=f"PAGA – pcs={n_pcs}, neighbours={n_neighbours}, counts={count_level}",
            show=False
        )

        plt.savefig(
            os.path.join(save_dir, f"paga_pcs{n_pcs}_neighbours{n_neighbours}_counts_{count_level}_celltypes.png"),
            bbox_inches="tight",
            dpi=150
        )
        plt.close()  # Clean up this figure immediately

        # Compute and plot FDG
        sc.tl.draw_graph(adata_test, init_pos="paga")
        sc.pl.embedding(
            adata_test,
            basis='draw_graph_fa',  # ← explicit basis (no X_ prefix needed)
            color="new_celltypes",
            title=f"FDG – pcs={n_pcs}, neighbours={n_neighbours}, counts={count_level}",
            legend_loc="right margin",
            legend_fontsize=8,
            size=18,
            frameon=False,  # optional but often nicer
            show=False
        )

        plt.savefig(
            os.path.join(save_dir, f"fdg_pcs{n_pcs}_neighbours{n_neighbours}_counts_{count_level}_celltypes.png"),
            bbox_inches="tight",
            dpi=150
        )
        plt.close()

        # write out the adata object for use in the RNA Velocity workflow
        adata_test.write_h5ad(
            os.path.join(save_dir, f"{sample_id}_fdg_pcs{n_pcs}_neighbours{n_neighbours}_counts_{count_level}_celltypes.h5ad"),
            compression="gzip"
        )
