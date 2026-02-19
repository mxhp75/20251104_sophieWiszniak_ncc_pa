#!/usr/bin/env python3
"""
SCVELO
=================

This script takes a filtered and normalised Seurat object and generates the intron/exon count files required for RNA-Velocity

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
from operator import index

#-----------------------------------------------------------------------------------
#   Global configs
#-----------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------
#   Import core libraries
#-----------------------------------------------------------------------------------

import pandas as pd
#import pandas_stubs
import scanpy as sc # single-cell analysis toolkit
import matplotlib.pyplot as plt
import os


#-----------------------------------------------------------------------------------
#   Params
#-----------------------------------------------------------------------------------
fPCA=50 # set no. of pcs input
fNEIGHBOUR=50 # set no. of neighbours input

# set the sample ID
sample_id = "e11_control"
base_dir = "/home/melanie-smith/workDir/sophieWiszniak/20251104_sophieWiszniak_ncc_pa/outDir/09-rna_velocity"
# set the sample id and output directory
out_dir = "/home/melanie-smith/workDir/sophieWiszniak/20251104_sophieWiszniak_ncc_pa/outDir/09-rna_velocity"

#-----------------------------------------------------------------------------------
# create embeddings for dynamic velocity modelling
# PAGA informed Force Directed Graph
#-----------------------------------------------------------------------------------

# add the adata file - this is the processed and filtered RNA count matrix from Seurat converted to .h5ad
adata = sc.read_h5ad(os.path.join(base_dir, sample_id, f"{sample_id}.h5ad"))
print(adata)
print(adata.layers.keys())       # at least 'counts', possibly others
print(adata.obsm.keys())         # X_pca, X_umap, etc.

# run pca
sc.tl.pca(adata, svd_solver='arpack', n_comps=fPCA)
# compute neighbourhood graph (constructs a k-nearest neighbour graph using PCA space)
sc.pp.neighbors(adata, n_neighbors=fNEIGHBOUR, n_pcs=fPCA)
# run paga new cell types (computes a coarse-grained connectivity graph between clusters)
sc.tl.paga(adata, groups='new_celltypes')
# Plot PAGA Graph
sc.pl.paga(adata, color=['new_celltypes'], node_size_scale=4, fontsize= 5) # plot paga graph

# create the PAGA informed FDG
sc.tl.draw_graph(adata, layout="fa", init_pos='paga')
sc.pl.draw_graph(adata, color="new_celltypes", legend_loc="right margin", legend_fontsize=8, size=18) # plot the FDG

pc_range = [10, 15, 20, 30]
neighbour_range = [15, 25, 30, 50]

# Create output directory if it doesn't exist
save_dir = os.path.join(out_dir, sample_id)
os.makedirs(save_dir, exist_ok=True)

for n_pcs in pc_range:
    for n_neighbours in neighbour_range:
        adata_test = adata.copy()

        sc.pp.neighbors(adata_test, n_pcs=n_pcs, n_neighbors=n_neighbours)
        sc.tl.paga(adata_test, groups="new_celltypes")
        sc.pl.paga(adata_test, threshold=0.03, show=False)
        sc.tl.draw_graph(adata_test, init_pos="paga")
        sc.pl.draw_graph(adata_test,
                         color="new_celltypes",
                         title=f"pcs={n_pcs}, neighbours={n_neighbours}",
                         show=False)

        plt.savefig(os.path.join(save_dir, f"fdg_pcs{n_pcs}_neighbours{n_neighbours}.png"),
                    bbox_inches="tight", dpi=150)
        plt.close()


# Save new paga informed FDG embedding and cell barcodes as a data frame
df_fa = pd.DataFrame(
    adata.obsm["X_draw_graph_fa"],
    columns=["FA_1", "FA_2"],
    index=adata.obs_names
).reset_index().rename(columns={"index": "barcode"})

df_fa.to_csv(os.path.join(save_dir, f"forceatlas2_embedding.csv"), index=False)