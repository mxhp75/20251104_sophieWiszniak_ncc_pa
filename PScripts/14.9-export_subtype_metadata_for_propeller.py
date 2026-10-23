#!/usr/bin/env python3
"""
Export per-cell metadata (sample, genotype, stage, endocardial subtype,
mesenchyme subtype) from the merged (all 4 samples) endocardium + mesenchyme
velocity object, for the R bootstrap propeller test on subtype proportions
(09-bootstrap-test-subtype-proportions.Rmd).

Endocardial subtype is read directly from 14.5's saved per-cell assignment
CSV. Mesenchyme subtype has no equivalent per-cell file (14.8 only saved the
genotype-level n/pct summary tables), so it is rebuilt here using the exact
same seurat_clusters -> name mapping as 14.8.

Author: Melanie Smith
Email:  melanie.smith@adelaide.edu.au
Date:   2026-10-21

Usage:
    python 14.9-export_subtype_metadata_for_propeller.py
"""

#-----------------------------------------------------------------------------------
#   Import core libraries
#-----------------------------------------------------------------------------------

import scanpy as sc
import pandas as pd
import numpy as np
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
velocity_dir = os.path.join(project_dir, "outDir", "14-rna_velocity")
egfp_dir = os.path.join(velocity_dir, "endocardial_mesenchymal_egfp")
out_dir = os.path.join(egfp_dir, "subtype_proportions_bootstrap")
os.makedirs(out_dir, exist_ok=True)

velocity_h5ad_path = os.path.join(
    velocity_dir, "velocity_output_dynamical", paga_group,
    f"all_samples_velocity_dynamical_{pcs}_{neighbours}_{count_level}_{paga_group}.h5ad"
)
endocardial_csv_path = os.path.join(egfp_dir, "all_samples_endocardial_subtype_assignments.csv")

# Sophie's seurat_clusters -> Mesenchyme sub-type mapping (same as 14.8)
cluster_to_name = {
    "0":  "EndoMT_derived",
    "2":  "EndoMT_derived (proliferating)",
    "3":  "E11_NCC",
    "8":  "E11_NCC (proliferating)",
    "5":  "E12_NCC_Vegfa",
    "10": "E12_NCC_Sema3c",
}

#-----------------------------------------------------------------------------------
#   Load merged velocity object and build the combined metadata table
#-----------------------------------------------------------------------------------

print("=== Loading merged velocity object ===")
adata = sc.read_h5ad(velocity_h5ad_path)

meta = adata.obs[["orig.ident", "genotype", "stage", "new_celltypes", "seurat_clusters"]].copy()

# ---- mesenchyme subtype (rebuilt, same logic as 14.8) ----
mes_mask = meta["new_celltypes"] == "Mesenchyme"
seurat_clusters_str = meta["seurat_clusters"].astype(str)
meta["mesenchyme_subtype"] = np.where(
    mes_mask, seurat_clusters_str.map(cluster_to_name).fillna("other"), "other"
)

# ---- endocardial subtype (read from 14.5's saved per-cell assignment) ----
endo_subtype = pd.read_csv(endocardial_csv_path, index_col=0)["endocardial_subtype"]
meta["endocardial_subtype"] = endo_subtype.reindex(meta.index)

#-----------------------------------------------------------------------------------
#   Sanity checks
#-----------------------------------------------------------------------------------

n_endo_missing = meta.loc[meta["new_celltypes"] == "Endocardium", "endocardial_subtype"].isna().sum()
print(f"\nEndocardium cells with missing subtype label: {n_endo_missing}")
assert n_endo_missing == 0, "Every Endocardium cell should have an endocardial_subtype label"

print("\nMesenchyme subtype counts:")
print(meta.loc[mes_mask, "mesenchyme_subtype"].value_counts())

print("\nEndocardial subtype counts:")
print(meta.loc[meta["new_celltypes"] == "Endocardium", "endocardial_subtype"].value_counts())

print("\nCell counts by genotype x stage:")
print(pd.crosstab(meta["genotype"], meta["stage"]))

#-----------------------------------------------------------------------------------
#   Save
#-----------------------------------------------------------------------------------

out_path = os.path.join(out_dir, "all_samples_subtype_metadata.csv")
meta.to_csv(out_path)
print(f"\nSaved: {out_path}")
