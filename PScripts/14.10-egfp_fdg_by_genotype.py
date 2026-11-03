#!/usr/bin/env python3
"""
EGFP+ cells on the Endocardium+Mesenchyme FDG, split by genotype - backup
sub-figure for the lab meeting deck (in case EGFP distribution comes up).

Uses the combined (all 4 samples merged) velocity object and its shared FDG
embedding - same embedding as the endocardial/mesenchyme subtype plots
(14.5/14.8), not separately-fitted per-genotype models - so the two genotype
panels are directly spatially comparable.

EGFP+ cells (618 of 14,631) are overwhelmingly in Mib1_KO (617) vs control
(1), as expected for a Cre-reporter linked to the KO allele.

Output: outDir/14-rna_velocity/endocardial_mesenchymal_egfp/egfp_enrichment/
  all_samples_egfp_fdg_by_genotype.png

Author: Melanie Smith
Email:  melanie.smith@adelaide.edu.au
Date:   2026-10-23

Usage:
    python 14.10-egfp_fdg_by_genotype.py
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

# EGFP_positive round-trips through h5ad as a Categorical (bool categories) -
# cast back to plain bool, then to a readable label for the legend
adata.obs["EGFP_positive"] = adata.obs["EGFP_positive"].astype(bool)
adata.obs["EGFP_status"] = np.where(adata.obs["EGFP_positive"], "EGFP+", "EGFP-")
adata.obs["EGFP_status"] = adata.obs["EGFP_status"].astype("category").cat.reorder_categories(
    ["EGFP-", "EGFP+"]
)

egfp_palette = ["#D3D3D3", "#2CA02C"]  # light grey, green - matches category order above

genotypes = ["control", "Mib1_KO"]

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for ax, geno in zip(axes, genotypes):
    sub = adata[adata.obs["genotype"] == geno]
    sc.pl.embedding(
        sub, basis="draw_graph_fa", color="EGFP_status", palette=egfp_palette,
        title=f"{geno} (n={sub.n_obs}, EGFP+={int(sub.obs['EGFP_positive'].sum())})",
        ax=ax, show=False, legend_loc="right margin" if geno == genotypes[-1] else None,
    )

fig.suptitle("EGFP+ cells on the Endocardium+Mesenchyme FDG, by genotype", fontsize=13)
fig.tight_layout()
fig.savefig(
    os.path.join(egfp_dir, "all_samples_egfp_fdg_by_genotype.png"),
    dpi=300, bbox_inches="tight"
)
plt.close(fig)

print("Saved: all_samples_egfp_fdg_by_genotype.png")
