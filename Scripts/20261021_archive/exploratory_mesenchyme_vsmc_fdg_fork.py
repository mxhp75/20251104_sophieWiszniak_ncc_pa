#!/usr/bin/env python3
"""
Exploratory fork - checking whether the Mesenchyme/VSMC gap seen on the
full-dataset UMAP (original and Harmony-integrated) is a real discontinuity
or a UMAP-specific exaggeration, via PAGA connectivity + a force-directed
graph (FDG) on the same neighbour graph.

Same recipe/parameters as 14.3-run_scanpy_paga_fdg_cell.py: neighbours graph
built on pcs=30/neighbours=25, PAGA grouped by cell type, FDG seeded from the
input embedding. Run twice - once on the original PCA embedding, once on the
Harmony embedding (both exported by 06.5c-export_embeddings_for_fdg_fork.Rmd)
- to see whether the two agree.

Standalone exploratory fork - reads only the CSVs exported by the R fork
script above, writes only to its own output folder. Does not touch any
existing script, object, or pipeline output.

Author: Melanie Smith
Email:  melanie.smith@adelaide.edu.au
Date:   2026-10-21

Usage:
    python exploratory_mesenchyme_vsmc_fdg_fork.py
"""

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

project_dir = "/home/melanie-smith/workDir/sophieWiszniak/20251104_sophieWiszniak_ncc_pa"
data_dir = os.path.join(
    project_dir, "outDir", "06.5-harmony_batch_integration_check", "mesenchyme_vsmc_fdg_check"
)
out_dir = data_dir

n_pcs = 30
n_neighbours = 25

meta = pd.read_csv(os.path.join(data_dir, "cell_metadata.csv"), index_col=0)
pca = pd.read_csv(os.path.join(data_dir, "pca_embeddings.csv"), index_col=0)
harmony = pd.read_csv(os.path.join(data_dir, "harmony_embeddings.csv"), index_col=0)

assert (meta.index == pca.index).all()
assert (meta.index == harmony.index).all()

adata = ad.AnnData(obs=meta)
adata.obsm["X_pca"] = pca.values
adata.obsm["X_harmony"] = harmony.values
adata.obs["new_celltypes"] = adata.obs["new_celltypes"].astype("category")

print(adata)

for rep_label, rep_key in [("original_pca", "X_pca"), ("harmony", "X_harmony")]:

    print(f"\n=== {rep_label} (pcs={n_pcs}, neighbours={n_neighbours}) ===")

    a = adata.copy()
    sc.pp.neighbors(a, n_pcs=n_pcs, n_neighbors=n_neighbours, use_rep=rep_key)
    sc.tl.paga(a, groups="new_celltypes")

    # PAGA connectivity matrix - the quantitative answer: how connected is
    # Mesenchyme to VSMC relative to other cell-type pairs
    conn = pd.DataFrame(
        a.uns["paga"]["connectivities"].toarray(),
        index=a.obs["new_celltypes"].cat.categories,
        columns=a.obs["new_celltypes"].cat.categories,
    )
    print("PAGA connectivity matrix:")
    print(conn.round(3))
    conn.to_csv(os.path.join(out_dir, f"paga_connectivity_{rep_label}.csv"))

    sc.pl.paga(
        a, threshold=0.01, color="new_celltypes", node_size_scale=4, fontsize=7,
        title=f"PAGA - {rep_label} - pcs={n_pcs}, neighbours={n_neighbours}", show=False,
    )
    plt.savefig(os.path.join(out_dir, f"paga_{rep_label}.png"), bbox_inches="tight", dpi=150)
    plt.close()

    sc.tl.draw_graph(a, init_pos=rep_key)
    sc.pl.embedding(
        a, basis="draw_graph_fa", color="new_celltypes",
        title=f"FDG - {rep_label} - pcs={n_pcs}, neighbours={n_neighbours}",
        legend_loc="right margin", legend_fontsize=8, size=10, frameon=False, show=False,
    )
    plt.savefig(os.path.join(out_dir, f"fdg_{rep_label}.png"), bbox_inches="tight", dpi=150)
    plt.close()

    a.write_h5ad(os.path.join(out_dir, f"fdg_{rep_label}.h5ad"), compression="gzip")

print("\nDone.")
