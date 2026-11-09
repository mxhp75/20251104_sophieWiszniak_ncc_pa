#!/usr/bin/env python3
"""
Follow-up plot requested at the 2026-10-26 meeting with Sophie, not added to
the additional_plots_meeting_20261013 lab book (saved locally only, to be
added by hand) - EGFP+ cells on the FDG, Mesenchyme cells only, split into 4
panels by orig.ident. Same green/light-grey convention as 14.10's
genotype-split version, just subset to Mesenchyme and split by sample
instead of genotype.

Uses the combined (all 4 samples merged) velocity object and its shared FDG
embedding - same embedding as 14.10/14.5/14.8 - so the four sample panels
are directly spatially comparable.

Output: outDir/14-rna_velocity/endocardial_mesenchymal_egfp/egfp_enrichment/
  all_samples_mesenchyme_egfp_fdg_by_sample.png

Author: Melanie Smith
Email:  melanie.smith@adelaide.edu.au
Date:   2026-11-02

Usage:
    python 14.17-mesenchyme_egfp_fdg_by_sample.py
"""

import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import os

project_dir = "/home/melanie-smith/workDir/sophieWiszniak/20251104_sophieWiszniak_ncc_pa"
out_dir = os.path.join(project_dir, "outDir", "14-rna_velocity")
egfp_dir = os.path.join(out_dir, "endocardial_mesenchymal_egfp", "egfp_enrichment")
os.makedirs(egfp_dir, exist_ok=True)

velocity_h5ad_path = os.path.join(
    out_dir, "velocity_output_dynamical", "celltypes",
    "all_samples_velocity_dynamical_pcs30_neighbours25_spliced_unspliced_celltypes.h5ad"
)

adata = sc.read_h5ad(velocity_h5ad_path)

adata.obs["EGFP_positive"] = adata.obs["EGFP_positive"].astype(bool)
adata.obs["EGFP_status"] = np.where(adata.obs["EGFP_positive"], "EGFP+", "EGFP-")
adata.obs["EGFP_status"] = adata.obs["EGFP_status"].astype("category").cat.reorder_categories(
    ["EGFP-", "EGFP+"]
)

egfp_palette = ["#D3D3D3", "#2CA02C"]  # light grey, green - matches category order above

mesenchyme = adata[adata.obs["new_celltypes"] == "Mesenchyme"]

samples = ["e11_control", "e11_ko", "e12_control", "e12_ko"]

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

for ax, sample in zip(axes.flat, samples):
    sub = mesenchyme[mesenchyme.obs["orig.ident"] == sample]
    sc.pl.embedding(
        sub, basis="draw_graph_fa", color="EGFP_status", palette=egfp_palette,
        title=f"{sample} (n={sub.n_obs}, EGFP+={int(sub.obs['EGFP_positive'].sum())})",
        ax=ax, show=False, legend_loc="right margin" if sample == samples[-1] else None,
    )

fig.suptitle("EGFP+ cells on the Mesenchyme-only FDG, by sample", fontsize=15)
fig.tight_layout()
fig.savefig(
    os.path.join(egfp_dir, "all_samples_mesenchyme_egfp_fdg_by_sample.png"),
    dpi=300, bbox_inches="tight"
)
plt.close(fig)

print("Saved: all_samples_mesenchyme_egfp_fdg_by_sample.png")
