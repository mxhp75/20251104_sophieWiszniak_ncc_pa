#!/usr/bin/env python3
"""
Endocardial subtype annotation (Sophie's marker scheme)
=================

Scores the four Endocardium subtype marker signatures Sophie described by
email (from her own inspection of the full merged dataset), checks them
against this project's genotype-level `13.x` clusters (`res_0.5`, chosen
specifically so the spliced/unspliced clustering resolves the same 4-way
split she is working with - see `13.2`), plots the raw signature scores,
then assigns each Endocardium cell to a subtype and re-plots by that label.

Signatures (Sophie's email, full-dataset cluster numbers in parentheses):
  - Ventricular endocardium (c15):        Irx5, Tgfbr3, Emcn, Frem1
  - Angiogenic/capillary endothelium (c16, E12.5-only):
                                           Apln, Col15a1, Kcne3, Clec1b,
                                           Fabp4, Dll4, Nrp2
  - Valve endocardium (c6):                Wnt9b, Adam23, Adamts8, Fam155a, Wnt4
  - Arterial/distal OFT endothelial (c11): Pak7, Fbln5, Rtl4, Isl1, Trpc5, Gja4

Verified against outDir/13-rna_velocity/{group}/{group}_spliced_endo_mesen_subset.h5ad
(the full-gene 13.2 object) before running this script: each signature scores
cleanly and near-exclusively in one res_0.5 cluster/cluster-pair, in both
genotypes - Ventricular in cluster 8, Angiogenic in cluster 9 (93-95% e12_*
cells, matching Sophie's "only appears at E12.5" note), Valve and Arterial
splitting the bulk cluster (control: 4 vs 5; ko: 5 vs 6 - cluster IDs are
arbitrary per independent clustering run, not directly comparable across
genotypes).

Two source objects per group, same pattern as `13.5`:
  - 13.2's full-gene .h5ad - scores the four signatures (full gene set needed)
  - 13.4's velocity .h5ad - target object with FDG/UMAP embeddings; scores
    merged on by cell barcode

Outputs per group (outDir/13-rna_velocity/{group}/endocardial_subtypes/):
  - signature_scores/  - FDG and UMAP FeaturePlots for each of the four raw
    signature scores
  - subtype_labels/    - FDG and UMAP plots coloured by the new
    `endocardial_subtype` column
  - {group}_endocardial_subtype_assignments.csv - per-cell barcode, the four
    scores, and the assigned subtype

New metadata column: `endocardial_subtype` - assigned only for Endocardium
cells (argmax of the four scores, or "Endocardium_unclassified" if all four
scores are <= 0); Mesenchyme cells keep their `new_celltypes` label
("Mesenchyme") unchanged. The existing `new_celltypes` column is left
untouched - this is an additional column, not a replacement.

Author: Melanie Smith
Email:  melanie.smith@adelaide.edu.au
Date:   2026-09-23

Usage:
    python 13.6-endocardial_subtype_annotation.py

Notes:
    - Main libraries used: scanpy, pandas, matplotlib
    - Python version: >= 3.9 recommended
"""

#-----------------------------------------------------------------------------------
#   Import core libraries
#-----------------------------------------------------------------------------------

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import warnings

#-----------------------------------------------------------------------------------
#   Params
#-----------------------------------------------------------------------------------

project_dir = "/home/melanie-smith/workDir/sophieWiszniak/20251104_sophieWiszniak_ncc_pa"
out_dir     = os.path.join(project_dir, "outDir", "13-rna_velocity")

group_names = ["control", "ko"]

paga_group  = "celltypes"
count_level = "spliced_unspliced"
pcs         = "pcs30"
neighbours  = "neighbours25"

signatures = {
    "Ventricular": ["Irx5", "Tgfbr3", "Emcn", "Frem1"],
    "Angiogenic":  ["Apln", "Col15a1", "Kcne3", "Clec1b", "Fabp4", "Dll4", "Nrp2"],
    "Valve":       ["Wnt9b", "Adam23", "Adamts8", "Fam155a", "Wnt4"],
    "Arterial":    ["Pak7", "Fbln5", "Rtl4", "Isl1", "Trpc5", "Gja4"],
}

#-----------------------------------------------------------------------------------
# Per-group processing
#-----------------------------------------------------------------------------------

assignment_summary_parts = []

for group_ID in group_names:

    print(f"\n=== Processing: {group_ID} ===")

    try:
        # ------------------------------------------------------------------
        # Score the four signatures on the full-gene 13.2 object
        # ------------------------------------------------------------------
        full_h5ad_path = os.path.join(out_dir, group_ID, f"{group_ID}_spliced_endo_mesen_subset.h5ad")
        adata_full = sc.read_h5ad(full_h5ad_path)

        for sig_name, genes in signatures.items():
            genes_present = [g for g in genes if g in adata_full.var_names]
            genes_missing = sorted(set(genes) - set(genes_present))
            if genes_missing:
                warnings.warn(f"{group_ID} — {sig_name}: missing genes {genes_missing}")
            sc.tl.score_genes(adata_full, gene_list=genes_present, score_name=sig_name)

        print(f"Scored {len(signatures)} signatures on {adata_full.n_vars} genes")

        # ------------------------------------------------------------------
        # Load the 13.4 velocity object and merge the signature scores onto it
        # ------------------------------------------------------------------
        velocity_h5ad_path = os.path.join(
            out_dir, group_ID, "velocity_output_dynamical", paga_group,
            f"{group_ID}_velocity_dynamical_{pcs}_{neighbours}_{count_level}_{paga_group}.h5ad"
        )
        adata = sc.read_h5ad(velocity_h5ad_path)

        for sig_name in signatures:
            adata.obs[sig_name] = adata_full.obs.loc[adata.obs_names, sig_name].values

        print(f"Merged signature scores onto velocity object ({adata.n_obs} cells)")

        # ------------------------------------------------------------------
        # Output directories
        # ------------------------------------------------------------------
        subtype_dir = os.path.join(out_dir, group_ID, "endocardial_subtypes")
        score_dir   = os.path.join(subtype_dir, "signature_scores")
        label_dir   = os.path.join(subtype_dir, "subtype_labels")
        for d in (score_dir, label_dir):
            os.makedirs(d, exist_ok=True)

        # ------------------------------------------------------------------
        # FDG and UMAP FeaturePlots of the four raw signature scores
        # ------------------------------------------------------------------
        for sig_name in signatures:
            for basis, basis_label in [("draw_graph_fa", "fdg"), ("umap_spliced", "umap")]:
                sc.pl.embedding(
                    adata, basis=basis, color=sig_name, cmap="viridis",
                    title=f"{group_ID} — {sig_name} — {basis_label}", show=False
                )
                plt.savefig(
                    os.path.join(score_dir, f"{group_ID}_{sig_name}_score_{basis_label}.png"),
                    dpi=300, bbox_inches="tight"
                )
                plt.close()

        print("Saved signature score FeaturePlots")

        # ------------------------------------------------------------------
        # Assign endocardial_subtype - argmax of the four scores for
        # Endocardium cells (or "Endocardium_unclassified" if all four scores
        # are <= 0), Mesenchyme cells keep their existing new_celltypes label.
        # new_celltypes itself is left untouched.
        # ------------------------------------------------------------------
        score_df = adata.obs[list(signatures.keys())]
        argmax_sig = score_df.idxmax(axis=1)
        max_score  = score_df.max(axis=1)

        endocardial_subtype = np.where(
            adata.obs["new_celltypes"] == "Endocardium",
            np.where(max_score > 0, argmax_sig, "Endocardium_unclassified"),
            adata.obs["new_celltypes"].astype(str)
        )
        adata.obs["endocardial_subtype"] = pd.Categorical(endocardial_subtype)

        counts = adata.obs["endocardial_subtype"].value_counts()
        print("endocardial_subtype counts:")
        print(counts)

        for subtype in ["Ventricular", "Angiogenic", "Valve", "Arterial", "Endocardium_unclassified"]:
            assignment_summary_parts.append({
                "genotype": group_ID,
                "endocardial_subtype": subtype,
                "n_cells": int(counts.get(subtype, 0)),
            })

        adata.obs[["latent_time"] + list(signatures.keys()) + ["endocardial_subtype"]].to_csv(
            os.path.join(subtype_dir, f"{group_ID}_endocardial_subtype_assignments.csv")
        )

        # ------------------------------------------------------------------
        # FDG and UMAP plots coloured by the new endocardial_subtype column
        # ------------------------------------------------------------------
        subtype_colours = {
            "Ventricular": "#984EA3",
            "Angiogenic":  "#FF7F00",
            "Valve":       "#377EB8",
            "Arterial":    "#E41A1C",
            "Endocardium_unclassified": "#999999",
            "Mesenchyme":  "#4DAF4A",
        }
        palette = [subtype_colours[c] for c in adata.obs["endocardial_subtype"].cat.categories]

        for basis, basis_label in [("draw_graph_fa", "fdg"), ("umap_spliced", "umap")]:
            sc.pl.embedding(
                adata, basis=basis, color="endocardial_subtype", palette=palette,
                title=f"{group_ID} — endocardial subtype — {basis_label}", show=False
            )
            plt.savefig(
                os.path.join(label_dir, f"{group_ID}_endocardial_subtype_{basis_label}.png"),
                dpi=300, bbox_inches="tight"
            )
            plt.close()

        print("Saved endocardial_subtype plots")

        print(f"Finished: {group_ID}")

    except Exception as e:
        plt.close("all")
        warnings.warn(f"{group_ID} — endocardial subtype annotation failed: {e}")

#-----------------------------------------------------------------------------------
# Combined subtype count summary across genotypes
#-----------------------------------------------------------------------------------

try:
    comparison_dir = os.path.join(out_dir, "bridge_markers_comparison", "endocardial_subtypes")
    os.makedirs(comparison_dir, exist_ok=True)

    summary_df = pd.DataFrame(assignment_summary_parts)
    summary_df.to_csv(
        os.path.join(comparison_dir, "endocardial_subtype_counts_by_genotype.csv"),
        index=False
    )
    print("\nEndocardial subtype counts by genotype:")
    print(summary_df.pivot(index="endocardial_subtype", columns="genotype", values="n_cells"))

except Exception as e:
    warnings.warn(f"Combined subtype summary failed: {e}")
