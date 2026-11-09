#!/usr/bin/env python3
"""
FDG with both Endocardium and Mesenchyme sub-clusters coloured together -
every cell in the Endocardium+Mesenchyme object gets its own subtype colour
(no grey "other" background), unlike 14.5/14.8's plots which each colour
one parent type's subtypes and grey out the other.

Reuses the already-exported per-cell subtype labels from 14.9
(all_samples_subtype_metadata.csv) and the same subtype colour palette as
14.5/14.8 (which matches the R-side propeller/bar-chart scripts).

Output: outDir/14-rna_velocity/endocardial_mesenchymal_egfp/subtype_labels/
  all_samples_combined_subtype_fdg.png / _umap.png

Author: Melanie Smith
Email:  melanie.smith@adelaide.edu.au
Date:   2026-10-26

Usage:
    python 14.12-combined_subtype_fdg.py
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

project_dir = "/home/melanie-smith/workDir/sophieWiszniak/20251104_sophieWiszniak_ncc_pa"
out_dir = os.path.join(project_dir, "outDir", "14-rna_velocity")
egfp_dir = os.path.join(out_dir, "endocardial_mesenchymal_egfp")
label_dir = os.path.join(egfp_dir, "subtype_labels")

velocity_h5ad_path = os.path.join(
    out_dir, "velocity_output_dynamical", "celltypes",
    "all_samples_velocity_dynamical_pcs30_neighbours25_spliced_unspliced_celltypes.h5ad"
)
subtype_csv_path = os.path.join(
    egfp_dir, "subtype_proportions_bootstrap", "all_samples_subtype_metadata.csv"
)

print("=== Loading velocity object and subtype metadata ===")
adata = sc.read_h5ad(velocity_h5ad_path)
subtype_meta = pd.read_csv(subtype_csv_path, index_col=0)

adata.obs["endocardial_subtype"] = subtype_meta["endocardial_subtype"].reindex(adata.obs_names)
adata.obs["mesenchyme_subtype"] = subtype_meta["mesenchyme_subtype"].reindex(adata.obs_names)

# every cell is Endocardium or Mesenchyme in this object - no "other" needed
is_endo = adata.obs["new_celltypes"] == "Endocardium"
combined_subtype = np.where(
    is_endo, adata.obs["endocardial_subtype"], adata.obs["mesenchyme_subtype"]
)
assert not pd.isna(combined_subtype).any(), "every cell should get a subtype label"

endocardial_order = ["Ventricular", "Angiogenic", "Valve", "Arterial"]
mesenchyme_order = [
    "EndoMT_derived", "EndoMT_derived (proliferating)",
    "E11_NCC", "E11_NCC (proliferating)",
    "E12_NCC_Vegfa", "E12_NCC_Sema3c",
]
combined_order = endocardial_order + mesenchyme_order

adata.obs["combined_subtype"] = pd.Categorical(combined_subtype, categories=combined_order)

# same palette as 14.5 (endocardial) + 14.8 (mesenchyme), which already
# matches the R-side propeller/bar-chart scripts
subtype_colours = {
    "Ventricular": "#66C2A5",
    "Angiogenic":  "#FC8D62",
    "Valve":       "#8DA0CB",
    "Arterial":    "#E78AC3",
    "EndoMT_derived":                 "#1B9E77",
    "EndoMT_derived (proliferating)": "#D95F02",
    "E11_NCC":                        "#7570B3",
    "E11_NCC (proliferating)":        "#E7298A",
    "E12_NCC_Vegfa":                  "#66A61E",
    "E12_NCC_Sema3c":                 "#E6AB02",
}
palette = [subtype_colours[c] for c in combined_order]

for basis, basis_label in [("draw_graph_fa", "fdg"), ("umap_spliced", "umap")]:
    ax = sc.pl.embedding(
        adata, basis=basis, color="combined_subtype", palette=palette,
        legend_loc=None, title="", show=False
    )
    legend_handles = [
        plt.Line2D([0], [0], marker="o", linestyle="", markersize=8, color=subtype_colours[c])
        for c in combined_order
    ]
    ax.legend(
        handles=legend_handles, labels=combined_order,
        loc="upper center", bbox_to_anchor=(0.5, -0.1),
        ncol=5, frameon=False, fontsize=9
    )
    plt.savefig(
        os.path.join(label_dir, f"all_samples_combined_subtype_{basis_label}.png"),
        dpi=300, bbox_inches="tight"
    )
    plt.close()

print("Saved combined subtype FDG/UMAP plots")
