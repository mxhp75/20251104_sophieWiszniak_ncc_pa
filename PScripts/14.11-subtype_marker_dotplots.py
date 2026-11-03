#!/usr/bin/env python3
"""
Subtype x marker-gene dot plots - compact evidence for the Endocardium and
Mesenchyme subtype calls, same format as the deck's existing slide 7
cell-type dot plot (dot size = % expressing, colour = average expression),
instead of showing the full grid of individual marker FeaturePlots from
15.7/15.8 (too many panels for a single slide).

Genes are Sophie's signature marker sets from 14.5 (signatures /
mesenchyme_signatures dicts). For the two "(proliferating)" Mesenchyme
subtypes, no separate identity gene set exists (14.5: "no distinct markers
given, so no separate signature") - shown instead against the same identity
genes as their non-proliferating counterpart, plus a Proliferation gene
block (Hist1h2ap; Esco2 dropped - not present in this object), so the plot
itself shows "same identity, elevated proliferation" rather than
duplicating identical gene columns under two different subtype labels.

Output: outDir/14-rna_velocity/endocardial_mesenchymal_egfp/
  subtype_labels/all_samples_endocardial_subtype_dotplot.png
  mesenchyme_subtypes_sophie/all_samples_mesenchyme_subtype_dotplot.png

Author: Melanie Smith
Email:  melanie.smith@adelaide.edu.au
Date:   2026-10-23

Usage:
    python 14.11-subtype_marker_dotplots.py
"""

import scanpy as sc
import pandas as pd
import numpy as np
import os

project_dir = "/home/melanie-smith/workDir/sophieWiszniak/20251104_sophieWiszniak_ncc_pa"
out_dir = os.path.join(project_dir, "outDir", "14-rna_velocity")
egfp_dir = os.path.join(out_dir, "endocardial_mesenchymal_egfp")

label_dir = os.path.join(egfp_dir, "subtype_labels")
mesenchyme_dir = os.path.join(egfp_dir, "mesenchyme_subtypes_sophie")

full_h5ad_path = os.path.join(out_dir, "all_samples_spliced_endo_mesen_subset.h5ad")
endocardial_csv_path = os.path.join(egfp_dir, "all_samples_endocardial_subtype_assignments.csv")

cluster_to_name = {
    "0":  "EndoMT_derived",
    "2":  "EndoMT_derived (proliferating)",
    "3":  "E11_NCC",
    "8":  "E11_NCC (proliferating)",
    "5":  "E12_NCC_Vegfa",
    "10": "E12_NCC_Sema3c",
}

# full gene set needed for a per-gene dot plot - the velocity object (used
# elsewhere for embeddings) has a reduced feature set and is missing some
# of these marker genes (e.g. Apln, Fam155a)
print("=== Loading full-gene-set object ===")
adata = sc.read_h5ad(full_h5ad_path)

# endocardial_subtype - read from 14.5's saved per-cell assignment
endo_subtype = pd.read_csv(endocardial_csv_path, index_col=0)["endocardial_subtype"]
adata.obs["endocardial_subtype"] = endo_subtype.reindex(adata.obs_names)

# mesenchyme_subtype - rebuilt, same logic as 14.8 (not saved to this h5ad)
mes_mask = adata.obs["new_celltypes"] == "Mesenchyme"
seurat_clusters_str = adata.obs["seurat_clusters"].astype(str)
adata.obs["mesenchyme_subtype"] = np.where(
    mes_mask, seurat_clusters_str.map(cluster_to_name).fillna("other"), "other"
)

#-----------------------------------------------------------------------------------
# Endocardial subtype dot plot
#-----------------------------------------------------------------------------------

endocardial_genes = {
    "Ventricular": ["Irx5", "Tgfbr3", "Emcn", "Frem1"],
    "Angiogenic":  ["Apln", "Col15a1", "Kcne3", "Clec1b", "Fabp4", "Dll4", "Nrp2"],
    "Valve":       ["Wnt9b", "Adam23", "Adamts8", "Fam155a", "Wnt4"],
    "Arterial":    ["Pak7", "Fbln5", "Rtl4", "Isl1", "Trpc5", "Gja4"],
}

import matplotlib.pyplot as plt

adata_endo = adata[adata.obs["new_celltypes"] == "Endocardium"].copy()
adata_endo.obs["endocardial_subtype"] = adata_endo.obs["endocardial_subtype"].astype("category")

sc.pl.dotplot(
    adata_endo, var_names=endocardial_genes, groupby="endocardial_subtype",
    standard_scale="var", cmap="Reds",
    show=False, save=False
)
plt.gcf().suptitle("Endocardial subtype marker genes", fontsize=16, y=1.12)
plt.savefig(
    os.path.join(label_dir, "all_samples_endocardial_subtype_dotplot.png"),
    dpi=300, bbox_inches="tight"
)
plt.close()
print("Saved endocardial subtype dot plot")

#-----------------------------------------------------------------------------------
# Mesenchyme subtype dot plot
#-----------------------------------------------------------------------------------

mesenchyme_genes = {
    "EndoMT_derived": ["Sall3", "Rapgef4", "Bmp6", "Htra1", "Lef1", "Nkain4"],
    "E11_NCC":        ["Plxna2", "Ramp1", "Kif26a", "Itga4"],
    "E12_NCC_Vegfa":  ["Vegfa", "Sema3d", "Ndufa4l2", "Chrm3", "Ctnna2"],
    "E12_NCC_Sema3c": ["Sema3c", "Cntn4", "Cntn6", "Syt1", "Penk"],
    "Proliferation":  ["Hist1h2ap"],
}

mesenchyme_order = [
    "EndoMT_derived", "EndoMT_derived (proliferating)",
    "E11_NCC", "E11_NCC (proliferating)",
    "E12_NCC_Vegfa", "E12_NCC_Sema3c",
]

adata_mes = adata[adata.obs["new_celltypes"] == "Mesenchyme"].copy()
adata_mes.obs["mesenchyme_subtype"] = pd.Categorical(
    adata_mes.obs["mesenchyme_subtype"], categories=mesenchyme_order
)

sc.pl.dotplot(
    adata_mes, var_names=mesenchyme_genes, groupby="mesenchyme_subtype",
    standard_scale="var", cmap="Reds",
    show=False, save=False
)
plt.gcf().suptitle("Mesenchyme subtype marker genes", fontsize=16, y=1.15)
plt.savefig(
    os.path.join(mesenchyme_dir, "all_samples_mesenchyme_subtype_dotplot.png"),
    dpi=300, bbox_inches="tight"
)
plt.close()
print("Saved mesenchyme subtype dot plot")
