#!/usr/bin/env python3
"""
Endocardial/mesenchymal subtype confirmation + EGFP enrichment (14.x)
=================

Narrower-scope follow-up to the `13.x` bridge marker analysis (`13.5`-`13.7`),
run on the `14.x` object (all 4 samples merged, cell-cycle genes retained,
no genotype split). The genotype-comparison parts of `13.x` don't apply here
(no genotype split), so this script covers only:

  1. Re-confirm Sophie's four Endocardium subtype signatures
     (Ventricular/Angiogenic/Valve/Arterial) hold up on this cc-retained,
     all-samples-merged object, independent of 13.x's no-cc/genotype-split
     design - cross-validation, not an assumption.
  2. Unsupervised DE characterisation of Mesenchyme res_0.4 sub-clusters
     (mirrors 13.7 - the existing marker panels don't separate Mesenchyme
     cleanly there either, not assumed to be true here without checking).
  3. NEW - not part of 13.x: does EGFP enrich in particular endocardial
     subtypes or particular Mesenchyme sub-clusters? Directly addresses
     whether EGFP marks a specific stage/subpopulation of the NCC-derived
     lineage, or Mesenchyme broadly.

Note: unlike 13.x (where res_0.5 cleanly resolved Valve from Arterial),
no resolution tested here (0.4-0.8) gives Valve its own distinct
sub-cluster - Arterial separates cleanly and consistently, but the largest
Endocardium cluster stays near-zero on the Valve score at every resolution
tested. The endocardial_subtype assignment below is per-cell (argmax of the
four continuous scores), so this doesn't block the analysis, but it is a
real, flagged difference from the 13.x result worth keeping in mind when
interpreting the Valve counts below.

Source objects:
  - 14.2's full-gene .h5ad (all_samples_spliced_endo_mesen_subset.h5ad) -
    scores the four Endocardium signatures (full gene set needed) and runs
    the Mesenchyme cluster DE (also needs the full gene set)
  - 14.4's velocity .h5ad - target object with FDG/UMAP embeddings and
    EGFP_UMI/EGFP_positive; scores merged on by cell barcode

Outputs (outDir/14-rna_velocity/endocardial_mesenchymal_egfp/):
  - signature_scores/      - FDG/UMAP FeaturePlots, four raw signature scores
  - subtype_labels/        - FDG/UMAP plots coloured by endocardial_subtype
  - mesenchyme_clusters/   - Mesenchyme res_0.4 cluster markers + highlight plots
  - egfp_enrichment/       - EGFP_UMI/EGFP_positive by endocardial_subtype and
    by Mesenchyme cluster (tables + plots) - the new analysis this script adds

Author: Melanie Smith
Email:  melanie.smith@adelaide.edu.au
Date:   2026-10-07

Usage:
    python 14.5-endocardial_mesenchymal_egfp_analysis.py

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
from matplotlib.colors import LinearSegmentedColormap
import os
import warnings

#-----------------------------------------------------------------------------------
#   Params
#-----------------------------------------------------------------------------------

project_dir = "/home/melanie-smith/workDir/sophieWiszniak/20251104_sophieWiszniak_ncc_pa"
out_dir     = os.path.join(project_dir, "outDir", "14-rna_velocity")
analysis_dir = os.path.join(out_dir, "endocardial_mesenchymal_egfp")

paga_group  = "celltypes"
count_level = "spliced_unspliced"
pcs         = "pcs30"
neighbours  = "neighbours25"

cluster_column = "res_0.4"  # 14.1/14.2's chosen resolution
min_cluster_size = 20       # exclude tiny stray Mesenchyme clusters from DE

signatures = {
    "Ventricular": ["Irx5", "Tgfbr3", "Emcn", "Frem1"],
    "Angiogenic":  ["Apln", "Col15a1", "Kcne3", "Clec1b", "Fabp4", "Dll4", "Nrp2"],
    "Valve":       ["Wnt9b", "Adam23", "Adamts8", "Fam155a", "Wnt4"],
    "Arterial":    ["Pak7", "Fbln5", "Rtl4", "Isl1", "Trpc5", "Gja4"],
}

# Sophie's proposed Mesenchyme sub-types (email, 2026-10-20), keyed to
# seurat_clusters cluster numbers (confirmed by direct marker-gene
# expression check - NOT the res_0.4 unsupervised clustering used in
# section 2 below, a different resolution/partition entirely):
#   cluster 0  -> EndoMT_derived
#   cluster 2  -> Proliferating (shares identity with cluster 0 + prolif markers)
#   cluster 3  -> E11_NCC
#   cluster 8  -> proliferating version of cluster 3 - no distinct markers given,
#                 so no separate signature (identifiable as E11_NCC + Proliferating
#                 both elevated)
#   cluster 5  -> E12_NCC_Vegfa
#   cluster 10 -> E12_NCC_Sema3c
mesenchyme_signatures = {
    "EndoMT_derived":  ["Sall3", "Rapgef4", "Bmp6", "Htra1", "Lef1", "Nkain4"],
    "Proliferating":   ["Hist1h2ap", "Esco2"],
    "E11_NCC":         ["Plxna2", "Ramp1", "Kif26a", "Itga4"],
    "E12_NCC_Vegfa":   ["Vegfa", "Sema3d", "Ndufa4l2", "Chrm3", "Ctnna2"],
    "E12_NCC_Sema3c":  ["Sema3c", "Cntn4", "Cntn6", "Syt1", "Penk"],
}

egfp_cmap = LinearSegmentedColormap.from_list("egfp", ["#e5e5e5", "gold", "red"])

#-----------------------------------------------------------------------------------
# Load objects, score signatures, merge onto the velocity object
#-----------------------------------------------------------------------------------

print("=== Loading objects ===")

full_h5ad_path = os.path.join(out_dir, "all_samples_spliced_endo_mesen_subset.h5ad")
adata_full = sc.read_h5ad(full_h5ad_path)

for sig_name, genes in signatures.items():
    genes_present = [g for g in genes if g in adata_full.var_names]
    genes_missing = sorted(set(genes) - set(genes_present))
    if genes_missing:
        warnings.warn(f"{sig_name}: missing genes {genes_missing}")
    sc.tl.score_genes(adata_full, gene_list=genes_present, score_name=sig_name)

print(f"Scored {len(signatures)} signatures on {adata_full.n_vars} genes")

for sig_name, genes in mesenchyme_signatures.items():
    genes_present = [g for g in genes if g in adata_full.var_names]
    genes_missing = sorted(set(genes) - set(genes_present))
    if genes_missing:
        warnings.warn(f"{sig_name}: missing genes {genes_missing}")
    sc.tl.score_genes(adata_full, gene_list=genes_present, score_name=sig_name)

print(f"Scored {len(mesenchyme_signatures)} Mesenchyme signatures on {adata_full.n_vars} genes")

velocity_h5ad_path = os.path.join(
    out_dir, "velocity_output_dynamical", paga_group,
    f"all_samples_velocity_dynamical_{pcs}_{neighbours}_{count_level}_{paga_group}.h5ad"
)
adata = sc.read_h5ad(velocity_h5ad_path)

# h5ad round-tripping stores EGFP_positive as a pandas Categorical (with bool
# categories) rather than a plain bool dtype - groupby().agg("sum") can't sum
# a Categorical directly, so cast it back explicitly.
adata.obs["EGFP_positive"] = adata.obs["EGFP_positive"].astype(bool)

for sig_name in signatures:
    adata.obs[sig_name] = adata_full.obs.loc[adata.obs_names, sig_name].values

for sig_name in mesenchyme_signatures:
    adata.obs[sig_name] = adata_full.obs.loc[adata.obs_names, sig_name].values

print(f"Merged signature scores onto velocity object ({adata.n_obs} cells)")

#-----------------------------------------------------------------------------------
# 1. Endocardial subtype confirmation
#-----------------------------------------------------------------------------------

try:
    score_dir = os.path.join(analysis_dir, "signature_scores")
    label_dir = os.path.join(analysis_dir, "subtype_labels")
    for d in (score_dir, label_dir):
        os.makedirs(d, exist_ok=True)

    for sig_name in signatures:
        for basis, basis_label in [("draw_graph_fa", "fdg"), ("umap_spliced", "umap")]:
            sc.pl.embedding(
                adata, basis=basis, color=sig_name, cmap="viridis",
                title=f"all_samples — {sig_name} — {basis_label}", show=False
            )
            plt.savefig(
                os.path.join(score_dir, f"all_samples_{sig_name}_score_{basis_label}.png"),
                dpi=300, bbox_inches="tight"
            )
            plt.close()

    print("Saved signature score FeaturePlots")

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

    adata.obs[["latent_time"] + list(signatures.keys()) + ["endocardial_subtype"]].to_csv(
        os.path.join(analysis_dir, "all_samples_endocardial_subtype_assignments.csv")
    )

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
            title=f"all_samples — endocardial subtype — {basis_label}", show=False
        )
        plt.savefig(
            os.path.join(label_dir, f"all_samples_endocardial_subtype_{basis_label}.png"),
            dpi=300, bbox_inches="tight"
        )
        plt.close()

    print("Saved endocardial_subtype plots")

except Exception as e:
    plt.close("all")
    warnings.warn(f"Endocardial subtype confirmation failed: {e}")

#-----------------------------------------------------------------------------------
# 1b. Mesenchyme sub-type signature scores (Sophie's seurat_clusters-based
# proposal, email 2026-10-20) - same signature_scores/ directory and
# all_samples_{signature}_score_{basis}.png naming as the Endocardium
# signatures above
#-----------------------------------------------------------------------------------

try:
    score_dir = os.path.join(analysis_dir, "signature_scores")
    os.makedirs(score_dir, exist_ok=True)

    for sig_name in mesenchyme_signatures:
        for basis, basis_label in [("draw_graph_fa", "fdg"), ("umap_spliced", "umap")]:
            sc.pl.embedding(
                adata, basis=basis, color=sig_name, cmap="viridis",
                title=f"all_samples — {sig_name} — {basis_label}", show=False
            )
            plt.savefig(
                os.path.join(score_dir, f"all_samples_{sig_name}_score_{basis_label}.png"),
                dpi=300, bbox_inches="tight"
            )
            plt.close()

    print("Saved Mesenchyme sub-type signature score plots")

except Exception as e:
    plt.close("all")
    warnings.warn(f"Mesenchyme sub-type signature scoring failed: {e}")

#-----------------------------------------------------------------------------------
# 2. Mesenchyme cluster characterisation (unsupervised, mirrors 13.7)
#-----------------------------------------------------------------------------------

try:
    mes_dir = os.path.join(analysis_dir, "mesenchyme_clusters")
    os.makedirs(mes_dir, exist_ok=True)

    mes_full = adata_full[adata_full.obs["new_celltypes"] == "Mesenchyme"].copy()

    cluster_sizes = mes_full.obs[cluster_column].value_counts()
    real_clusters = cluster_sizes[cluster_sizes >= min_cluster_size].index.tolist()
    excluded = cluster_sizes[cluster_sizes < min_cluster_size]
    print(f"\n{cluster_column} Mesenchyme clusters: {cluster_sizes.to_dict()}")
    if len(excluded):
        print(f"Excluding tiny stray clusters (< {min_cluster_size} cells): {excluded.to_dict()}")

    mes_real = mes_full[mes_full.obs[cluster_column].isin(real_clusters)].copy()
    mes_real.obs[cluster_column] = mes_real.obs[cluster_column].cat.remove_unused_categories()

    sc.tl.rank_genes_groups(
        mes_real, groupby=cluster_column, groups=real_clusters,
        method="wilcoxon", pts=True
    )

    marker_rows = []
    for cl in real_clusters:
        de = sc.get.rank_genes_groups_df(mes_real, group=cl)
        sig = de[(de["pvals_adj"] < 0.05) & (de["pct_nz_group"] > 0.25)]
        sig_up = sig[sig["logfoldchanges"] > 0].sort_values("scores", ascending=False)
        top = sig_up.head(15).copy()
        top.insert(0, "cluster", cl)
        marker_rows.append(top)
        print(f"\ncluster {cl} ({cluster_sizes[cl]} cells) — top 10 markers:")
        print(top.head(10)[["names", "logfoldchanges", "pvals_adj", "pct_nz_group"]].to_string(index=False))

    marker_df = pd.concat(marker_rows, ignore_index=True)
    marker_df.to_csv(os.path.join(mes_dir, "all_samples_mesenchyme_cluster_markers.csv"), index=False)

    sc.pl.rank_genes_groups_dotplot(
        mes_real, groupby=cluster_column, n_genes=5, show=False,
        title=f"all_samples — Mesenchyme {cluster_column} cluster markers"
    )
    plt.savefig(
        os.path.join(mes_dir, "all_samples_mesenchyme_cluster_marker_dotplot.png"),
        dpi=300, bbox_inches="tight"
    )
    plt.close()

    adata.obs["mesenchyme_cluster_plot"] = np.where(
        adata.obs["new_celltypes"] == "Mesenchyme",
        adata.obs[cluster_column].astype(str), "other"
    )
    adata.obs["mesenchyme_cluster_plot"] = np.where(
        adata.obs["mesenchyme_cluster_plot"].isin(real_clusters),
        adata.obs["mesenchyme_cluster_plot"], "other"
    )
    adata.obs["mesenchyme_cluster_plot"] = pd.Categorical(adata.obs["mesenchyme_cluster_plot"])

    cmap = plt.get_cmap("tab10")
    mes_palette = {cl: cmap(i % 10) for i, cl in enumerate(real_clusters)}
    mes_palette["other"] = (0.85, 0.85, 0.85, 1.0)
    colours = [mes_palette[c] for c in adata.obs["mesenchyme_cluster_plot"].cat.categories]

    for basis, basis_label in [("draw_graph_fa", "fdg"), ("umap_spliced", "umap")]:
        sc.pl.embedding(
            adata, basis=basis, color="mesenchyme_cluster_plot", palette=colours,
            title=f"all_samples — Mesenchyme {cluster_column} clusters — {basis_label}", show=False
        )
        plt.savefig(
            os.path.join(mes_dir, f"all_samples_mesenchyme_clusters_{basis_label}.png"),
            dpi=300, bbox_inches="tight"
        )
        plt.close()

    print("Saved Mesenchyme cluster characterisation outputs")

except Exception as e:
    plt.close("all")
    warnings.warn(f"Mesenchyme cluster characterisation failed: {e}")

#-----------------------------------------------------------------------------------
# 3. EGFP enrichment by endocardial subtype and by Mesenchyme cluster - new
#-----------------------------------------------------------------------------------

try:
    egfp_dir = os.path.join(analysis_dir, "egfp_enrichment")
    os.makedirs(egfp_dir, exist_ok=True)

    # --- by endocardial subtype (includes Mesenchyme as its own category) ---
    by_subtype = adata.obs.groupby("endocardial_subtype").agg(
        n_cells=("EGFP_UMI", "size"),
        n_egfp_positive=("EGFP_positive", "sum"),
        mean_EGFP_UMI=("EGFP_UMI", "mean"),
    )
    by_subtype["pct_egfp_positive"] = 100 * by_subtype["n_egfp_positive"] / by_subtype["n_cells"]
    by_subtype = by_subtype.sort_values("pct_egfp_positive", ascending=False)
    by_subtype.to_csv(os.path.join(egfp_dir, "egfp_enrichment_by_endocardial_subtype.csv"))
    print("\nEGFP enrichment by endocardial_subtype:")
    print(by_subtype)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(by_subtype.index.astype(str), by_subtype["pct_egfp_positive"],
           color=[subtype_colours.get(c, "grey") for c in by_subtype.index])
    ax.set_ylabel("% EGFP positive")
    ax.set_title("EGFP positivity by endocardial subtype")
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
    fig.tight_layout()
    fig.savefig(os.path.join(egfp_dir, "egfp_pct_positive_by_endocardial_subtype.png"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)

    # --- by Mesenchyme cluster specifically ---
    mes_mask = (adata.obs["new_celltypes"] == "Mesenchyme").values
    mes_cluster_series = adata.obs.loc[mes_mask, cluster_column].astype(str)
    mes_cluster_series = np.where(mes_cluster_series.isin(real_clusters), mes_cluster_series, "other")

    mes_egfp_df = pd.DataFrame({
        "mesenchyme_cluster": mes_cluster_series,
        "EGFP_UMI": adata.obs.loc[mes_mask, "EGFP_UMI"].values,
        "EGFP_positive": adata.obs.loc[mes_mask, "EGFP_positive"].values,
    })
    by_mes_cluster = mes_egfp_df.groupby("mesenchyme_cluster").agg(
        n_cells=("EGFP_UMI", "size"),
        n_egfp_positive=("EGFP_positive", "sum"),
        mean_EGFP_UMI=("EGFP_UMI", "mean"),
    )
    by_mes_cluster["pct_egfp_positive"] = 100 * by_mes_cluster["n_egfp_positive"] / by_mes_cluster["n_cells"]
    by_mes_cluster = by_mes_cluster.sort_values("pct_egfp_positive", ascending=False)
    by_mes_cluster.to_csv(os.path.join(egfp_dir, "egfp_enrichment_by_mesenchyme_cluster.csv"))
    print("\nEGFP enrichment by Mesenchyme cluster:")
    print(by_mes_cluster)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(by_mes_cluster.index.astype(str), by_mes_cluster["pct_egfp_positive"], color="#4DAF4A", alpha=0.8)
    ax.set_xlabel(f"Mesenchyme {cluster_column} cluster")
    ax.set_ylabel("% EGFP positive")
    ax.set_title("EGFP positivity by Mesenchyme cluster")
    fig.tight_layout()
    fig.savefig(os.path.join(egfp_dir, "egfp_pct_positive_by_mesenchyme_cluster.png"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)

    print("\nSaved EGFP enrichment outputs")

except Exception as e:
    plt.close("all")
    warnings.warn(f"EGFP enrichment analysis failed: {e}")

print("\nFinished: 14.5 endocardial/mesenchymal/EGFP analysis")
