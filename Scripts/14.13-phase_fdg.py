#!/usr/bin/env python3
"""
FDG of the Endocardium+Mesenchyme merged object (all 4 samples), coloured
by cell cycle phase (G1/S/G2M) rather than subtype.

Uses the same phase colour convention established elsewhere in this
project (e.g. the deck's "Cell cycle is a strong signal" slide and
08-bootstrap-test-celltype-proportions-cell-cycle.Rmd's phase_cols) -
G1 = blue, S = orange, G2M = green - not the subtype palette, since this
is a different categorical variable.

Output: outDir/14-rna_velocity/endocardial_mesenchymal_egfp/
  phase/all_samples_phase_fdg.png / _umap.png

Author: Melanie Smith
Email:  melanie.smith@adelaide.edu.au
Date:   2026-10-26

Usage:
    python 14.13-phase_fdg.py
"""

import scanpy as sc
import matplotlib.pyplot as plt
import os

project_dir = "/home/melanie-smith/workDir/sophieWiszniak/20251104_sophieWiszniak_ncc_pa"
out_dir = os.path.join(project_dir, "outDir", "14-rna_velocity")
egfp_dir = os.path.join(out_dir, "endocardial_mesenchymal_egfp")
phase_dir = os.path.join(egfp_dir, "phase")
os.makedirs(phase_dir, exist_ok=True)

velocity_h5ad_path = os.path.join(
    out_dir, "velocity_output_dynamical", "celltypes",
    "all_samples_velocity_dynamical_pcs30_neighbours25_spliced_unspliced_celltypes.h5ad"
)

print("=== Loading merged velocity object ===")
adata = sc.read_h5ad(velocity_h5ad_path)

# same phase colours used throughout this project (R's phase_cols:
# G1 = "#1f77b4", S = "#ff7f0e", G2M = "#2ca02c")
phase_order = ["G1", "S", "G2M"]
phase_colours = {
    "G1":  "#1f77b4",
    "S":   "#ff7f0e",
    "G2M": "#2ca02c",
}
adata.obs["Phase"] = adata.obs["Phase"].astype("category").cat.reorder_categories(phase_order)
palette = [phase_colours[c] for c in phase_order]

for basis, basis_label in [("draw_graph_fa", "fdg"), ("umap_spliced", "umap")]:
    ax = sc.pl.embedding(
        adata, basis=basis, color="Phase", palette=palette,
        legend_loc=None, title="", show=False
    )
    legend_handles = [
        plt.Line2D([0], [0], marker="o", linestyle="", markersize=8, color=phase_colours[c])
        for c in phase_order
    ]
    ax.legend(
        handles=legend_handles, labels=phase_order,
        loc="upper center", bbox_to_anchor=(0.5, -0.1),
        ncol=3, frameon=False, fontsize=9
    )
    plt.savefig(
        os.path.join(phase_dir, f"all_samples_phase_{basis_label}.png"),
        dpi=300, bbox_inches="tight"
    )
    plt.close()

print("Saved cell cycle phase FDG/UMAP plots")
