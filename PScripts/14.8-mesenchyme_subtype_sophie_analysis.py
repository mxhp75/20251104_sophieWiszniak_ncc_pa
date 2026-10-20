#!/usr/bin/env python3
"""
Mesenchyme sub-types proposed by Sophie, confirmed against seurat_clusters
(email, 2026-10-20) - formalises the ad hoc analysis run for the lab
meeting deck follow-up.

Sophie sent six Mesenchyme sub-cluster numbers with proposed identities
and marker genes. Confirmed these are seurat_clusters cluster IDs (not the
res_0.4 unsupervised clustering 14.5 already uses for the Mesenchyme
cluster/EGFP analysis - a different partition of the same cells; see the
14.1 lab book, section 6b, for the full identification writeup) by
checking that every marker set peaks in its claimed cluster. Two of the
six are proliferating counterparts of two others (cluster 2 proliferates
cluster 0's identity, cluster 8 proliferates cluster 3's), so there are
four distinct identities plus their cycling states:

  cluster 0  -> EndoMT_derived
  cluster 2  -> EndoMT_derived (proliferating)
  cluster 3  -> E11_NCC
  cluster 8  -> E11_NCC (proliferating)
  cluster 5  -> E12_NCC_Vegfa
  cluster 10 -> E12_NCC_Sema3c

Not a new dynamical-model fit - reads the already-computed 14.4 velocity
.h5ad (seurat_clusters, genotype and new_celltypes are already columns on
it) and only adds the categorical mesenchyme_subtype label plus a
%-by-genotype breakdown. The five continuous signature scores
(EndoMT_derived/Proliferating/E11_NCC/E12_NCC_Vegfa/E12_NCC_Sema3c) that
these category names echo were already added to 14.5
(signature_scores/all_samples_{signature}_score_{fdg,umap}.png) - this
script is the categorical, per-cell-label counterpart, plus the
genotype breakdown 14.5 doesn't compute.

Outputs (outDir/14-rna_velocity/endocardial_mesenchymal_egfp/mesenchyme_subtypes_sophie/):
  - all_samples_mesenchyme_subtype_sophie_fdg.png / _umap.png - categorical
    label on both bases, full object (other cell types shown grey)
  - mesenchyme_subtype_n_by_genotype.csv / mesenchyme_subtype_pct_by_genotype.csv -
    counts and % of Mesenchyme cells in each sub-type, control vs KO.
    Computed directly from the single shared seurat_clusters column (not
    from separately-fitted per-genotype models, unlike the 13.x
    endocardial-subtype comparison) - the clustering itself doesn't
    change between genotypes here, so this is arguably more robust than
    comparing across two independently-fit models. No statistical test
    is run - propeller bootstrap would be the natural next step, same as
    flagged for the endocardial subtypes.

Author: Melanie Smith
Email:  melanie.smith@adelaide.edu.au
Date:   2026-10-20

Usage:
    python 14.8-mesenchyme_subtype_sophie_analysis.py
"""

#-----------------------------------------------------------------------------------
#   Import core libraries
#-----------------------------------------------------------------------------------

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

#-----------------------------------------------------------------------------------
#   Params
#-----------------------------------------------------------------------------------

paga_group = "celltypes"
count_level = "spliced_unspliced"
pcs_value = 30
neighbours_value = 25
pcs = f"pcs{pcs_value}"
neighbours = f"neighbours{neighbours_value}"

project_dir = "/home/melanie-smith/workDir/sophieWiszniak/20251104_sophieWiszniak_ncc_pa"
out_dir     = os.path.join(project_dir, "outDir", "14-rna_velocity")
analysis_dir = os.path.join(out_dir, "endocardial_mesenchymal_egfp", "mesenchyme_subtypes_sophie")
os.makedirs(analysis_dir, exist_ok=True)

velocity_h5ad_path = os.path.join(
    out_dir, "velocity_output_dynamical", paga_group,
    f"all_samples_velocity_dynamical_{pcs}_{neighbours}_{count_level}_{paga_group}.h5ad"
)

# Sophie's seurat_clusters -> Mesenchyme sub-type mapping (see docstring)
cluster_to_name = {
    "0":  "EndoMT_derived",
    "2":  "EndoMT_derived (proliferating)",
    "3":  "E11_NCC",
    "8":  "E11_NCC (proliferating)",
    "5":  "E12_NCC_Vegfa",
    "10": "E12_NCC_Sema3c",
}

subtype_colours = {
    "EndoMT_derived":                "#E41A1C",
    "EndoMT_derived (proliferating)": "#FF7F00",
    "E11_NCC":                        "#4DAF4A",
    "E11_NCC (proliferating)":        "#984EA3",
    "E12_NCC_Vegfa":                  "#377EB8",
    "E12_NCC_Sema3c":                 "#A65628",
    "other":                          (0.85, 0.85, 0.85, 1.0),
}

#-----------------------------------------------------------------------------------
# Build the categorical mesenchyme_subtype label
#-----------------------------------------------------------------------------------

print("=== Loading velocity object ===")
adata = sc.read_h5ad(velocity_h5ad_path)

mes_mask = adata.obs["new_celltypes"] == "Mesenchyme"
seurat_clusters_str = adata.obs["seurat_clusters"].astype(str)
subtype = np.where(mes_mask, seurat_clusters_str.map(cluster_to_name).fillna("other"), "other")
adata.obs["mesenchyme_subtype"] = pd.Categorical(subtype)

print("mesenchyme_subtype counts:")
print(adata.obs.loc[mes_mask, "mesenchyme_subtype"].value_counts())

#-----------------------------------------------------------------------------------
# % of Mesenchyme cells in each sub-type, by genotype
#-----------------------------------------------------------------------------------

n_by_genotype = pd.crosstab(
    adata.obs.loc[mes_mask, "mesenchyme_subtype"], adata.obs.loc[mes_mask, "genotype"]
)
pct_by_genotype = 100 * n_by_genotype / n_by_genotype.sum(axis=0)
pct_by_genotype = pct_by_genotype.sort_values("control", ascending=False)

print("\n% of Mesenchyme cells by sub-type and genotype:")
print(pct_by_genotype.round(1))

n_by_genotype.to_csv(os.path.join(analysis_dir, "mesenchyme_subtype_n_by_genotype.csv"))
pct_by_genotype.to_csv(os.path.join(analysis_dir, "mesenchyme_subtype_pct_by_genotype.csv"))

#-----------------------------------------------------------------------------------
# Categorical FDG / UMAP plots, full object (other cell types shown grey) -
# same style as 14.5's endocardial_subtype plot
#-----------------------------------------------------------------------------------

palette = [subtype_colours[c] for c in adata.obs["mesenchyme_subtype"].cat.categories]

for basis, basis_label in [("draw_graph_fa", "fdg"), ("umap_spliced", "umap")]:
    sc.pl.embedding(
        adata, basis=basis, color="mesenchyme_subtype", palette=palette,
        title=f"all_samples - Mesenchyme subtype (Sophie) - {basis_label}", show=False
    )
    plt.savefig(
        os.path.join(analysis_dir, f"all_samples_mesenchyme_subtype_sophie_{basis_label}.png"),
        dpi=300, bbox_inches="tight"
    )
    plt.close()

print("\nSaved mesenchyme_subtype plots and genotype breakdown CSVs")
print("Finished: 14.8 Mesenchyme subtype (Sophie) analysis")
