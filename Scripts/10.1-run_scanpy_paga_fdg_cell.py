#!/usr/bin/env python3
"""

=================

This script takes a filtered and normalised Seurat object, converted previously to a .h5ad file, and generates the paga graph (coarse-grained connectivity map between cell type clusters) and uses this information generate the FDG projections

Author: Melanie Smith
Email:  melanie.smith@adelaide.edu.au
Date:   2026-09-02

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

# set the sample IDs
sample_ids = [
    "e11_control",
    "e12_control",
    "e11_ko",
    "e12_ko"
]

base_dir = "/home/melanie-smith/workDir/sophieWiszniak/20251104_sophieWiszniak_ncc_pa/outDir/10-rna_velocity"
# set the sample id and output directory
out_dir = "/home/melanie-smith/workDir/sophieWiszniak/20251104_sophieWiszniak_ncc_pa/outDir/10-rna_velocity"

# set the count level (ie gene or spliced/unspliced)
#count_level = "gene"
count_level = "spliced_unspliced"

pc_range = [10, 15, 20, 30] # this is how many PCs will be considered -> slices X_pca_no_cc_sample (50 dims available)
neighbour_range = [15, 25, 30, 50]

#-----------------------------------------------------------------------------------
# Sanity check — plot the Seurat-computed UMAP for each sample
# Compare these against your Seurat UMAPs before proceeding to PAGA/FDG
#-----------------------------------------------------------------------------------

print("\nRunning UMAP sanity checks...\n")

for sample_id in sample_ids:

    print(f"\n{'='*60}")
    print(f"Sanity check - sample: {sample_id}")
    print(f"{'='*60}\n")

    adata = sc.read_h5ad(os.path.join(out_dir, sample_id, f"{sample_id}_{count_level}.h5ad"))

    sc.pl.embedding(
        adata,
        basis='X_umap_no_cc_sample',
        color='new_celltypes',
        title=f"{sample_id} - umap_no_cc_sample",
        legend_loc='right margin',
        show=False
    )

    plt.savefig(
        os.path.join(out_dir, sample_id, f"{sample_id}_umap_no_cc_sample_check.png"),
        bbox_inches='tight',
        dpi=150
    )
    plt.close()

    print(f"Plot saved to: {out_dir}/{sample_id}/{sample_id}_umap_no_cc_sample_check.png")

print("\nSanity checks complete. Please review the UMAP plots before proceeding.\n")

#-----------------------------------------------------------------------------------
# Loop over all samples
# PAGA informed Force Directed Graph
#-----------------------------------------------------------------------------------

for sample_id in sample_ids:

    print(f"\n{'='*60}")
    print(f"Processing sample: {sample_id}")
    print(f"{'='*60}\n")

    #-----------------------------------------------------------------------------------
    # create embeddings for dynamic velocity modelling
    # PAGA informed Force Directed Graph
    #-----------------------------------------------------------------------------------

    # add the adata file - this is the processed and filtered RNA count matrix from Seurat converted to .h5ad
    # it contains gene level and spliced/unspliced level UMAP embeddings
    adata = sc.read_h5ad(os.path.join(out_dir, sample_id, f"{sample_id}_{count_level}.h5ad"))

    # let's take a look at the object
    print(adata)
    print(adata.layers.keys())       # at least 'counts', possibly others
    print(adata.obsm.keys())         # X_pca, X_umap, X_umap_spliced

    # REMOVED: sc.tl.pca() call removed to preserve the pre-computed cell-cycle-regressed
    # PCA from Seurat (X_pca_no_cc_sample). Running sc.tl.pca() here would overwrite
    # adata.obsm['X_pca'] with a fresh scanpy PCA computed from raw counts, discarding
    # the CC gene removal and variable feature selection performed in Seurat.

    #-----------------------------------------------------------------------------------
    # Explore the nPC and nNeighbours range
    # PAGA informed Force Directed Graph
    # Save output for all permutations
    #-----------------------------------------------------------------------------------

    # Create output directory if it doesn't exist
    save_dir = os.path.join(out_dir, sample_id, f"scanpy_output/celltypes")
    os.makedirs(save_dir, exist_ok=True)

    for n_pcs in pc_range:
        for n_neighbours in neighbour_range:
            adata_test = adata.copy()

            # use_rep points to the Seurat-computed, CC-regressed, sample-level PCA
            # n_pcs slices the first n dimensions of X_pca_no_cc_sample (50 dims available)
            # select the appropriate PCA embedding based on count level
            use_rep = 'X_pca_no_cc_sample' if count_level == "gene" else 'X_pca_spliced'

            sc.pp.neighbors(adata_test, n_pcs=n_pcs, n_neighbors=n_neighbours, use_rep=use_rep)

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

    print(f"\nFinished sample: {sample_id}\n")

print("\nAll samples complete.\n")