#!/usr/bin/env python3
"""
Bridge marker analysis
=================

Investigates the sparse "bridge" of cells connecting the Endocardium and
Mesenchyme populations in the 13.x (cell-cycle-genes-removed) RNA velocity
result, using four curated marker panels:
  - Bridge_Endocardial:     Ptprb, Egfl7, Bmx, Icam2, Irx6, Cytl1, Ecscr,
                             Cdh5, Tie1, Rasip1, Plxnd1
  - Bridge_NCMesenchyme:    Sox10, Twist1, Ednra, Sema3c, Prrx1, Prrx2, Sox9
  - Bridge_EndMTMesenchyme: Snai1, Snai2, Twist2, Cdh11, Has2, Tbx20
  - Bridge_Fibroblast:      Tcf21, Pdgfra, Col1a1, Col3a1, Dcn, Lum, Postn, Fn1

Rather than gating "bridge" cells by eye on FDG/UMAP coordinates (not
reproducible across re-runs), cells are characterised along the already-
computed `latent_time` axis (root = Endocardium, t=0).

Two source objects are used per group:
  - 13.2's full-gene .h5ad (~24.6k genes) - used ONLY to score the four
    marker panels, since three genes (Irx6, Cytl1, Tcf21) are absent from
    13.4's final velocity object, which is subset to the top 2000 dispersion
    genes used for the dynamical model fit.
  - 13.4's velocity .h5ad (2000 genes) - the target object, carrying
    `latent_time`, the FDG/UMAP embeddings, and the existing 13.1 EndMT_*
    module scores. The four new panel scores are merged onto this object by
    cell barcode (confirmed identical barcode sets between 13.2 and 13.4).

Outputs per group (outDir/13-rna_velocity/{group}/bridge_markers/):
  - latent_time_trends/  - mean +/- SEM of each panel score across ~20
    latent-time quantile bins, plus a summary CSV
  - featureplots_fdg/    - FDG FeaturePlots for the four panel scores and
    four key individual markers (Cdh5, Ecscr, Twist1, Postn)
  - latent_time_diagnostics/ - follow-up on the ko latent-time trend bending
    back upward past ~0.7: flags cells with latent_time > 0.7 and
    Bridge_Endocardial > 0.3 ("reversal" cells), reports which cluster(s)
    (see `cluster_column`) they belong to, and highlights them on FDG/UMAP
  - extremity_cluster_characterisation/ - follow-up on the reversal cells:
    identifies the cluster(s) accounting for the bulk of them (a
    coherent subpopulation sitting at a topological extremity/spur of the
    Endocardium body in both genotypes, not scattered noise), then runs a
    Wilcoxon DE test (full cluster membership vs rest of Endocardium, on
    the full 13.2 gene set) plus a QC metric comparison, a marker dotplot,
    and FDG/UMAP highlight plots of the full cluster membership, plus a
    direct comparison of the four Bridge_* panel scores (extremity cluster
    vs rest of Endocardium) - confirms the extremity cluster is MORE
    endocardial and LESS mesenchymal/EndMT-like than average, i.e. it is
    not a transitional/bridge population despite its elevated latent_time

Combined outputs (outDir/13-rna_velocity/bridge_markers_comparison/):
  - zone_comparison/     - cells split into latent-time terciles
    (Endocardial-committed / Bridge / Mesenchymal-committed), panel scores
    compared across zone and genotype
  - latent_time_diagnostics/ - reversal-cluster summary and latent_time
    distribution comparison across genotypes

Author: Melanie Smith
Email:  melanie.smith@adelaide.edu.au
Date:   2026-09-23

Usage:
    python 13.5-bridge_marker_analysis.py

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

n_latent_time_bins = 20

# The 13.2 spliced/unspliced-level clustering column used to identify Endocardium
# sub-clusters (e.g. the extremity cluster below). Must match 13.2's chosen_resolutions
# ("res_0.5" for both groups, revised from "res_0.4" - at 0.4 the spliced/unspliced
# clustering lumps valve and arterial endocardium together in both genotypes, whereas
# 0.5 cleanly splits them into the 4-way endocardial scheme Sophie is working with).
# Named as a single constant rather than hardcoded per-use, since cluster IDs are
# arbitrary per independent clustering run and this project has already been bitten
# once by a stale hardcoded resolution string.
cluster_column = "res_0.5"

marker_panels = {
    "Bridge_Endocardial":     ["Ptprb", "Egfl7", "Bmx", "Icam2", "Irx6", "Cytl1",
                                "Ecscr", "Cdh5", "Tie1", "Rasip1", "Plxnd1"],
    "Bridge_NCMesenchyme":    ["Sox10", "Twist1", "Ednra", "Sema3c", "Prrx1", "Prrx2", "Sox9"],
    "Bridge_EndMTMesenchyme": ["Snai1", "Snai2", "Twist2", "Cdh11", "Has2", "Tbx20"],
    "Bridge_Fibroblast":      ["Tcf21", "Pdgfra", "Col1a1", "Col3a1", "Dcn", "Lum", "Postn", "Fn1"],
}

# Existing EndMT module scores from 13.1, carried through to 13.4 unchanged -
# included in the trend plot for continuity with the earlier analysis
existing_endmt_scores = ["EndMT_endocardial1", "EndMT_transition1",
                          "EndMT_mesenchyme1", "EndMT_Notch1", "EndMT_TGFb1"]

key_individual_markers = ["Cdh5", "Ecscr", "Twist1", "Postn"]

zone_labels = ["Endocardial-committed", "Bridge", "Mesenchymal-committed"]

# Thresholds for flagging "reversal" cells - high latent_time despite a high
# Bridge_Endocardial score - chosen from visual inspection of the latent-time
# trend plots (where the ko aggregate trend bends back upward past ~0.7)
reversal_lt_threshold    = 0.7
reversal_score_threshold = 0.3

#-----------------------------------------------------------------------------------
# Per-group processing
#-----------------------------------------------------------------------------------

combined_zone_df_parts     = []
reversal_summary_parts     = []
extremity_size_summary_parts = []

for group_ID in group_names:

    print(f"\n=== Processing: {group_ID} ===")

    try:
        # ------------------------------------------------------------------
        # Score the four marker panels on the full-gene 13.2 object
        # ------------------------------------------------------------------
        full_h5ad_path = os.path.join(out_dir, group_ID, f"{group_ID}_spliced_endo_mesen_subset.h5ad")
        adata_full = sc.read_h5ad(full_h5ad_path)

        for panel_name, genes in marker_panels.items():
            genes_present = [g for g in genes if g in adata_full.var_names]
            genes_missing = sorted(set(genes) - set(genes_present))
            if genes_missing:
                warnings.warn(f"{group_ID} — {panel_name}: missing genes {genes_missing}")
            sc.tl.score_genes(adata_full, gene_list=genes_present, score_name=panel_name)

        print(f"Scored {len(marker_panels)} marker panels on {adata_full.n_vars} genes")

        # ------------------------------------------------------------------
        # Load the 13.4 velocity object and merge the panel scores onto it
        # ------------------------------------------------------------------
        velocity_h5ad_path = os.path.join(
            out_dir, group_ID, "velocity_output_dynamical", paga_group,
            f"{group_ID}_velocity_dynamical_{pcs}_{neighbours}_{count_level}_{paga_group}.h5ad"
        )
        adata = sc.read_h5ad(velocity_h5ad_path)

        for panel_name in marker_panels:
            adata.obs[panel_name] = adata_full.obs.loc[adata.obs_names, panel_name].values

        print(f"Merged panel scores onto velocity object ({adata.n_obs} cells)")

        # ------------------------------------------------------------------
        # Set up output directories
        # ------------------------------------------------------------------
        bridge_dir  = os.path.join(out_dir, group_ID, "bridge_markers")
        trend_dir   = os.path.join(bridge_dir, "latent_time_trends")
        feature_dir = os.path.join(bridge_dir, "featureplots_fdg")
        for d in (trend_dir, feature_dir):
            os.makedirs(d, exist_ok=True)

        # ------------------------------------------------------------------
        # Latent time trend plot - mean +/- SEM of panel scores across bins
        # ------------------------------------------------------------------
        score_cols = list(marker_panels.keys()) + existing_endmt_scores
        df = adata.obs[["latent_time"] + score_cols].copy()
        df["lt_bin"] = pd.qcut(df["latent_time"], q=n_latent_time_bins, duplicates="drop")
        df["lt_bin_mid"] = df["lt_bin"].apply(lambda iv: iv.mid).astype(float)

        summary = df.groupby("lt_bin_mid")[score_cols].agg(["mean", "sem"])
        summary.columns = ["_".join(col).strip("_") for col in summary.columns]
        summary = summary.reset_index()

        summary.to_csv(
            os.path.join(bridge_dir, f"{group_ID}_bridge_marker_latent_time_summary.csv"),
            index=False
        )

        # Full trend plot: four new panels + five existing EndMT scores
        fig, ax = plt.subplots(figsize=(11, 6))
        colors = plt.cm.tab10(np.linspace(0, 1, len(score_cols)))
        for i, col in enumerate(score_cols):
            ax.plot(summary["lt_bin_mid"], summary[f"{col}_mean"], label=col, color=colors[i])
            ax.fill_between(
                summary["lt_bin_mid"],
                summary[f"{col}_mean"] - summary[f"{col}_sem"],
                summary[f"{col}_mean"] + summary[f"{col}_sem"],
                alpha=0.15, color=colors[i]
            )
        ax.set_xlabel("Latent time")
        ax.set_ylabel("Mean module score")
        ax.set_title(f"{group_ID} — bridge marker panels vs. latent time")
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
        fig.tight_layout()
        fig.savefig(
            os.path.join(trend_dir, f"{group_ID}_bridge_marker_latent_time_trend.png"),
            dpi=300, bbox_inches="tight"
        )
        plt.close(fig)

        # Cleaner version: the four new panels only
        fig, ax = plt.subplots(figsize=(9, 6))
        colors4 = plt.cm.tab10(np.linspace(0, 1, len(marker_panels)))
        for i, col in enumerate(marker_panels):
            ax.plot(summary["lt_bin_mid"], summary[f"{col}_mean"], label=col, color=colors4[i], linewidth=2)
            ax.fill_between(
                summary["lt_bin_mid"],
                summary[f"{col}_mean"] - summary[f"{col}_sem"],
                summary[f"{col}_mean"] + summary[f"{col}_sem"],
                alpha=0.15, color=colors4[i]
            )
        ax.set_xlabel("Latent time")
        ax.set_ylabel("Mean module score")
        ax.set_title(f"{group_ID} — bridge marker panels (new panels only) vs. latent time")
        ax.legend(fontsize=9)
        fig.tight_layout()
        fig.savefig(
            os.path.join(trend_dir, f"{group_ID}_bridge_marker_panels_only_latent_time_trend.png"),
            dpi=300, bbox_inches="tight"
        )
        plt.close(fig)

        print("Saved latent time trend plots")

        # ------------------------------------------------------------------
        # FDG FeaturePlots - four panel scores + four key individual markers
        # ------------------------------------------------------------------
        for panel_name in marker_panels:
            sc.pl.embedding(
                adata, basis="draw_graph_fa", color=panel_name, cmap="viridis",
                title=f"{group_ID} — {panel_name} — FDG", show=False
            )
            plt.savefig(
                os.path.join(feature_dir, f"{group_ID}_{panel_name}_fdg.png"),
                dpi=300, bbox_inches="tight"
            )
            plt.close()

        for gene in key_individual_markers:
            if gene not in adata.var_names:
                warnings.warn(f"{group_ID} — {gene} not in the 2000-gene velocity object, skipping FDG FeaturePlot")
                continue
            sc.pl.embedding(
                adata, basis="draw_graph_fa", color=gene, cmap="viridis",
                title=f"{group_ID} — {gene} — FDG", show=False
            )
            plt.savefig(
                os.path.join(feature_dir, f"{group_ID}_{gene}_fdg.png"),
                dpi=300, bbox_inches="tight"
            )
            plt.close()

        print("Saved FDG FeaturePlots")

        # ------------------------------------------------------------------
        # Assign latent-time tercile zones, stash for the combined
        # cross-genotype comparison after the group loop
        # ------------------------------------------------------------------
        zone_df = adata.obs[["latent_time"] + list(marker_panels.keys())].copy()
        zone_df["zone"] = pd.qcut(zone_df["latent_time"], q=3, labels=zone_labels)
        zone_df["genotype"] = group_ID
        combined_zone_df_parts.append(zone_df)

        # ------------------------------------------------------------------
        # Latent-time diagnostics - do late-latent-time, high-Bridge_Endocardial
        # ("reversal") cells map to a genuine spatial/cluster subpopulation,
        # or are they scattered noise? Follow-up investigation prompted by
        # the ko latent-time trend bending back upward past ~0.7.
        # ------------------------------------------------------------------
        lt_diag_dir = os.path.join(bridge_dir, "latent_time_diagnostics")
        os.makedirs(lt_diag_dir, exist_ok=True)

        reversal_mask = (
            (adata.obs["latent_time"] > reversal_lt_threshold) &
            (adata.obs["Bridge_Endocardial"] > reversal_score_threshold)
        )
        adata.obs["reversal_flag"] = np.where(reversal_mask, "reversal", "other")

        n_late     = int((adata.obs["latent_time"] > reversal_lt_threshold).sum())
        n_reversal = int(reversal_mask.sum())
        pct_of_late = 100 * n_reversal / n_late if n_late > 0 else np.nan

        print(
            f"Latent-time diagnostics: {n_reversal} reversal cells "
            f"(latent_time>{reversal_lt_threshold} & Bridge_Endocardial>{reversal_score_threshold}), "
            f"of {n_late} cells with latent_time>{reversal_lt_threshold} total "
            f"({pct_of_late:.1f}% are reversal cells)"
        )

        reversal_summary_parts.append({
            "genotype": group_ID,
            "n_cells": adata.n_obs,
            "n_late_latent_time": n_late,
            "n_reversal_cells": n_reversal,
            "pct_of_late_that_are_reversal": pct_of_late,
        })

        cluster_breakdown = adata.obs.loc[reversal_mask, cluster_column].value_counts()
        cluster_breakdown.to_csv(os.path.join(lt_diag_dir, f"{group_ID}_reversal_cell_cluster_breakdown.csv"))

        for basis, basis_label in [("draw_graph_fa", "fdg"), ("umap_spliced", "umap")]:
            coords_key = f"X_{basis}" if f"X_{basis}" in adata.obsm else basis
            coords = adata.obsm[coords_key]

            fig, ax = plt.subplots(figsize=(8, 6))
            ax.scatter(coords[~reversal_mask.values, 0], coords[~reversal_mask.values, 1],
                       s=6, c="lightgray", label="other")
            ax.scatter(coords[reversal_mask.values, 0], coords[reversal_mask.values, 1],
                       s=10, c="red",
                       label=f"reversal (lt>{reversal_lt_threshold} & BridgeEndo>{reversal_score_threshold})")
            ax.set_title(f"{group_ID} — reversal cells highlighted on {basis_label}")
            ax.legend(fontsize=8)
            ax.axis("off")
            fig.savefig(
                os.path.join(lt_diag_dir, f"{group_ID}_reversal_highlight_{basis_label}.png"),
                dpi=300, bbox_inches="tight"
            )
            plt.close(fig)

        print("Saved latent-time diagnostic outputs")

        # ------------------------------------------------------------------
        # Extremity-cluster characterisation - what makes the graph-extremity
        # Endocardium subcluster (identified above via the reversal cells)
        # different from the rest of Endocardium? The cluster(s) are
        # identified dynamically as whichever cluster(s) (in `cluster_column`)
        # account for the large majority of reversal cells, rather than
        # hardcoding cluster labels (arbitrary per-group, since clustering is
        # run independently for control and ko, and can also shift between
        # resolutions - e.g. this population was clusters 7/8 at res_0.4 and
        # is clusters 8/9 at res_0.5).
        # ------------------------------------------------------------------
        extremity_dir = os.path.join(bridge_dir, "extremity_cluster_characterisation")
        os.makedirs(extremity_dir, exist_ok=True)

        cluster_frac = (cluster_breakdown / cluster_breakdown.sum()).sort_values(ascending=False)
        n_needed = int((cluster_frac.cumsum() < 0.9).sum()) + 1
        extremity_clusters = cluster_frac.index[:n_needed].tolist()

        # Explicit clustering-column label - cluster IDs are arbitrary per
        # independent clustering run and can shift between resolutions, so
        # always state which column they refer to. Do not confuse with 13.1's
        # "seurat_clusters_endo_mesen" gene-level clustering, which can label
        # this same cell population differently (confirmed via cell-level
        # overlap when this was res_0.4 vs seurat_clusters_endo_mesen).
        extremity_cluster_label = f"{cluster_column} clusters {', '.join(extremity_clusters)}"

        print(
            f"Extremity cluster(s) for {group_ID}: {extremity_cluster_label} "
            f"(covering {100 * cluster_frac.iloc[:n_needed].sum():.1f}% of reversal cells)"
        )

        endo_full = adata_full[adata_full.obs["new_celltypes"] == "Endocardium"].copy()
        endo_full.obs["extremity_flag"] = pd.Categorical(np.where(
            endo_full.obs[cluster_column].isin(extremity_clusters), "extremity_cluster", "rest_of_endocardium"
        ))

        n_extremity = int((endo_full.obs["extremity_flag"] == "extremity_cluster").sum())
        n_rest = int((endo_full.obs["extremity_flag"] == "rest_of_endocardium").sum())
        print(f"  {n_extremity} extremity-cluster cells vs {n_rest} rest-of-Endocardium cells")

        extremity_size_summary_parts.append({
            "genotype": group_ID,
            "n_extremity_cluster": n_extremity,
            "n_rest_of_endocardium": n_rest,
            "n_total_endocardium": n_extremity + n_rest,
            "pct_of_endocardium_in_extremity_cluster": 100 * n_extremity / (n_extremity + n_rest),
        })

        # QC metric comparison (rule out a trivial low-quality/doublet explanation)
        qc_cols = [c for c in ["nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo",
                                "scDblFinder.score", "scDblFinder.weighted"]
                   if c in endo_full.obs.columns]
        qc_summary = endo_full.obs.groupby("extremity_flag")[qc_cols].mean().T
        qc_summary.to_csv(os.path.join(extremity_dir, f"{group_ID}_extremity_cluster_qc_comparison.csv"))

        # Differential expression: extremity cluster vs rest of Endocardium
        sc.tl.rank_genes_groups(
            endo_full, groupby="extremity_flag", groups=["extremity_cluster"],
            reference="rest_of_endocardium", method="wilcoxon", pts=True
        )
        de = sc.get.rank_genes_groups_df(endo_full, group="extremity_cluster")
        sig = de[(de["pvals_adj"] < 0.05) & (de["pct_nz_group"] > 0.25)].copy()
        sig_up = sig[sig["logfoldchanges"] > 0].sort_values("scores", ascending=False)
        sig_down = sig[sig["logfoldchanges"] < 0].sort_values("scores", ascending=True)

        sig_up.to_csv(os.path.join(extremity_dir, f"{group_ID}_extremity_cluster_UP_markers.csv"), index=False)
        sig_down.to_csv(os.path.join(extremity_dir, f"{group_ID}_extremity_cluster_DOWN_markers.csv"), index=False)
        print(
            f"  {len(sig_up)} significant up genes, {len(sig_down)} significant down genes "
            f"(padj<0.05, expressed in >25% of extremity-cluster cells)"
        )

        sc.pl.rank_genes_groups_dotplot(
            endo_full, groupby="extremity_flag", n_genes=12, show=False,
            title=f"{group_ID} — {extremity_cluster_label} vs rest of Endocardium"
        )
        plt.savefig(
            os.path.join(extremity_dir, f"{group_ID}_extremity_cluster_marker_dotplot.png"),
            dpi=300, bbox_inches="tight"
        )
        plt.close()

        # Full cluster membership (not just the reversal-flagged subset) highlighted
        # on the velocity object's FDG/UMAP embeddings
        extremity_barcodes = set(endo_full.obs_names[endo_full.obs["extremity_flag"] == "extremity_cluster"])
        extremity_membership_mask = adata.obs_names.isin(extremity_barcodes)

        for basis, basis_label in [("draw_graph_fa", "fdg"), ("umap_spliced", "umap")]:
            coords_key = f"X_{basis}" if f"X_{basis}" in adata.obsm else basis
            coords = adata.obsm[coords_key]

            fig, ax = plt.subplots(figsize=(8, 6))
            ax.scatter(coords[~extremity_membership_mask, 0], coords[~extremity_membership_mask, 1],
                       s=6, c="lightgray", label="other")
            ax.scatter(coords[extremity_membership_mask, 0], coords[extremity_membership_mask, 1],
                       s=10, c="darkorange", label=f"{extremity_cluster_label} (full membership)")
            ax.set_title(f"{group_ID} — {extremity_cluster_label} on {basis_label}")
            ax.legend(fontsize=8)
            ax.axis("off")
            fig.savefig(
                os.path.join(extremity_dir, f"{group_ID}_extremity_cluster_highlight_{basis_label}.png"),
                dpi=300, bbox_inches="tight"
            )
            plt.close(fig)

        print("Saved extremity-cluster characterisation outputs")

        # ------------------------------------------------------------------
        # Does the extremity cluster actually look like a bridge/transitional
        # population? Compare its own Bridge_* panel scores against the rest
        # of Endocardium (not just latent_time) - this is the direct test of
        # whether the elevated latent_time reflects real progression toward
        # Mesenchyme, or a cluster that is if anything MORE purely endocardial
        # and LESS mesenchymal than average.
        # ------------------------------------------------------------------
        endo_velocity_mask = (adata.obs["new_celltypes"] == "Endocardium").values
        panel_cols = list(marker_panels.keys())

        bridge_score_df = adata.obs.loc[endo_velocity_mask, panel_cols].copy()
        bridge_score_df["group"] = np.where(
            extremity_membership_mask[endo_velocity_mask], "extremity_cluster", "rest_of_endocardium"
        )
        bridge_score_comparison = bridge_score_df.groupby("group")[panel_cols].mean()
        bridge_score_comparison.to_csv(
            os.path.join(extremity_dir, f"{group_ID}_extremity_cluster_bridge_panel_scores.csv")
        )
        print("Bridge panel scores, extremity cluster vs rest of Endocardium:")
        print(bridge_score_comparison)

        fig, ax = plt.subplots(figsize=(8, 5))
        bridge_score_comparison.T.plot(kind="bar", ax=ax, color=["darkorange", "gray"])
        ax.axhline(0, color="black", linewidth=0.8)
        ax.set_ylabel("Mean module score")
        ax.set_title(f"{group_ID} — {extremity_cluster_label}: bridge panel scores vs rest of Endocardium")
        ax.legend(title="")
        plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
        fig.tight_layout()
        fig.savefig(
            os.path.join(extremity_dir, f"{group_ID}_extremity_cluster_bridge_panel_scores.png"),
            dpi=300, bbox_inches="tight"
        )
        plt.close(fig)

        print("Saved extremity-cluster bridge-panel-score comparison")

        print(f"Finished: {group_ID}")

    except Exception as e:
        plt.close("all")
        warnings.warn(f"{group_ID} — bridge marker analysis failed: {e}")

#-----------------------------------------------------------------------------------
# Combined cross-genotype zone comparison
#-----------------------------------------------------------------------------------

try:
    comparison_dir = os.path.join(out_dir, "bridge_markers_comparison", "zone_comparison")
    os.makedirs(comparison_dir, exist_ok=True)

    combined_zone_df = pd.concat(combined_zone_df_parts, ignore_index=True)

    combined_zone_df.to_csv(
        os.path.join(comparison_dir, "bridge_marker_zone_scores_by_genotype.csv"),
        index=False
    )

    # Fraction of cells per zone per genotype
    zone_fractions = (
        combined_zone_df.groupby("genotype")["zone"]
        .value_counts(normalize=True)
        .rename("fraction")
        .reset_index()
    )
    zone_fractions.to_csv(
        os.path.join(comparison_dir, "bridge_zone_fractions_by_genotype.csv"),
        index=False
    )
    print("\nZone fractions by genotype:")
    print(zone_fractions)

    # Violin plots: each panel score across zone, split by genotype
    for panel_name in marker_panels:

        fig, axes = plt.subplots(1, 1, figsize=(10, 5.5))
        positions = []
        zone_centers = []
        data_to_plot = []
        tick_labels = []
        colors_by_genotype = {"control": "#377EB8", "ko": "#E41A1C"}

        pos = 0
        for zone in zone_labels:
            zone_start = pos
            for genotype in group_names:
                subset = combined_zone_df[
                    (combined_zone_df["zone"] == zone) & (combined_zone_df["genotype"] == genotype)
                ][panel_name]
                data_to_plot.append(subset.values)
                positions.append(pos)
                tick_labels.append(genotype)
                pos += 1
            zone_centers.append((zone_start + pos - 1) / 2)
            pos += 1.2

        parts = axes.violinplot(data_to_plot, positions=positions, showmeans=True, showextrema=False)
        for i, body in enumerate(parts["bodies"]):
            genotype = group_names[i % len(group_names)]
            body.set_facecolor(colors_by_genotype[genotype])
            body.set_alpha(0.7)

        axes.set_xticks(positions)
        axes.set_xticklabels(tick_labels, fontsize=9)

        # Zone name labels, one per zone group, above the genotype ticks
        y_min, y_max = axes.get_ylim()
        label_y = y_min - 0.14 * (y_max - y_min)
        for zone, center in zip(zone_labels, zone_centers):
            axes.text(center, label_y, zone, ha="center", va="top", fontsize=10, fontweight="bold")
        axes.set_ylim(y_min - 0.22 * (y_max - y_min), y_max)

        axes.set_ylabel(panel_name)
        axes.set_title(f"{panel_name} — by latent-time zone and genotype")
        fig.tight_layout()
        fig.savefig(
            os.path.join(comparison_dir, f"{panel_name}_zone_comparison_violin.png"),
            dpi=300, bbox_inches="tight"
        )
        plt.close(fig)

    print(f"Completed zone comparison for {len(marker_panels)} panels")

except Exception as e:
    plt.close("all")
    warnings.warn(f"Combined zone comparison failed: {e}")

#-----------------------------------------------------------------------------------
# Combined latent-time diagnostics - reversal-cluster summary + distribution
# comparison across genotypes (follow-up to the per-group reversal highlight
# plots above)
#-----------------------------------------------------------------------------------

try:
    lt_diag_comparison_dir = os.path.join(out_dir, "bridge_markers_comparison", "latent_time_diagnostics")
    os.makedirs(lt_diag_comparison_dir, exist_ok=True)

    reversal_summary_df = pd.DataFrame(reversal_summary_parts)
    reversal_summary_df.to_csv(
        os.path.join(lt_diag_comparison_dir, "reversal_cluster_summary_by_genotype.csv"),
        index=False
    )
    print("\nReversal-cluster summary by genotype:")
    print(reversal_summary_df)

    lt_summary = combined_zone_df.groupby("genotype")["latent_time"].describe()
    lt_summary.to_csv(os.path.join(lt_diag_comparison_dir, "latent_time_distribution_summary_by_genotype.csv"))
    print("\nLatent time distribution by genotype:")
    print(lt_summary)

    colors_by_genotype = {"control": "#377EB8", "ko": "#E41A1C"}

    fig, ax = plt.subplots(figsize=(8, 5))
    for genotype in group_names:
        vals = combined_zone_df.loc[combined_zone_df["genotype"] == genotype, "latent_time"]
        ax.hist(vals, bins=30, density=True, alpha=0.5, color=colors_by_genotype[genotype], label=genotype)
    ax.set_xlabel("Latent time")
    ax.set_ylabel("Density")
    ax.set_title("Latent time distribution by genotype")
    ax.legend()
    fig.tight_layout()
    fig.savefig(
        os.path.join(lt_diag_comparison_dir, "latent_time_distribution_by_genotype.png"),
        dpi=300, bbox_inches="tight"
    )
    plt.close(fig)

    print("Saved latent time distribution comparison")

except Exception as e:
    plt.close("all")
    warnings.warn(f"Combined latent-time diagnostics failed: {e}")

#-----------------------------------------------------------------------------------
# Combined extremity-cluster size comparison - is the graph-extremity
# Endocardium subpopulation itself a larger fraction of Endocardium in ko,
# independent of the latent-time-tail dilution effect above?
#-----------------------------------------------------------------------------------

try:
    extremity_comparison_dir = os.path.join(out_dir, "bridge_markers_comparison", "extremity_cluster_characterisation")
    os.makedirs(extremity_comparison_dir, exist_ok=True)

    extremity_size_df = pd.DataFrame(extremity_size_summary_parts)
    extremity_size_df.to_csv(
        os.path.join(extremity_comparison_dir, "extremity_cluster_size_by_genotype.csv"),
        index=False
    )
    print("\nExtremity cluster size by genotype:")
    print(extremity_size_df)

    colors_by_genotype = {"control": "#377EB8", "ko": "#E41A1C"}

    fig, ax = plt.subplots(figsize=(6, 5))
    bar_colors = [colors_by_genotype[g] for g in extremity_size_df["genotype"]]
    ax.bar(extremity_size_df["genotype"], extremity_size_df["pct_of_endocardium_in_extremity_cluster"],
           color=bar_colors, alpha=0.8)
    for i, row in extremity_size_df.iterrows():
        ax.text(i, row["pct_of_endocardium_in_extremity_cluster"] + 0.3,
                f"{row['n_extremity_cluster']}/{row['n_total_endocardium']}",
                ha="center", fontsize=9)
    ax.set_ylabel("% of Endocardium in extremity cluster")
    ax.set_title("Extremity-cluster size as a fraction of Endocardium, by genotype")
    fig.tight_layout()
    fig.savefig(
        os.path.join(extremity_comparison_dir, "extremity_cluster_size_by_genotype.png"),
        dpi=300, bbox_inches="tight"
    )
    plt.close(fig)

    print("Saved extremity cluster size comparison")

except Exception as e:
    plt.close("all")
    warnings.warn(f"Combined extremity-cluster size comparison failed: {e}")
