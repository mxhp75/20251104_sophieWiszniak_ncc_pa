#!/usr/bin/env python3
"""
Mesenchymal cluster characterisation (unsupervised)
=================

Follow-up to `13.6`'s Endocardium subtype annotation. The three existing
Bridge_* mesenchymal marker panels (from `13.5`) were checked against
Mesenchyme's `res_0.5` sub-clusters before writing this script, and do not
cleanly separate them: `Bridge_Fibroblast` dominates every single
sub-cluster in both genotypes (pairwise score correlations all low-positive,
0.14-0.24), unlike Endocardium where each cluster was dominated by a
different one of Sophie's four signatures. An argmax-style labeling would
therefore just call nearly every Mesenchyme cell "Fibroblast" - not a
meaningful subtype split.

This script instead takes an unsupervised approach: for each genotype's
existing `res_0.5` Mesenchyme sub-clusters (already computed in `13.2`,
excluding tiny strays with < 20 cells), run a one-vs-rest Wilcoxon DE test
to discover what actually distinguishes each cluster, rather than assuming
which marker axis (origin / transition-state / fibroblast-type) matters.
This mirrors `13.5`'s extremity-cluster-characterisation approach (find
markers, then interpret), not `13.6`'s marker-panel-argmax approach.

No subtype names are assigned in this script - that needs domain
interpretation of the resulting marker genes (the same way Sophie's email
was needed to name the Endocardium clusters), not a guess.

Source object per group: 13.2's full-gene .h5ad (needed for DE across all
~24.6k genes, not just the 2000-gene 13.4 velocity object).

Outputs per group (outDir/13-rna_velocity/{group}/mesenchymal_subtypes/):
  - cluster_characterisation/{group}_mesenchyme_cluster_markers.csv - top
    significant up-regulated genes per cluster (one-vs-rest, padj<0.05,
    expressed in >25% of the cluster)
  - cluster_characterisation/{group}_mesenchyme_cluster_marker_dotplot.png -
    summary dotplot, top 5 genes per cluster across all clusters
  - cluster_characterisation/{group}_mesenchyme_clusters_{fdg,umap}.png -
    Mesenchyme coloured by res_0.5 cluster, Endocardium shown in grey for
    spatial context
  - arterial_connectivity/{group}_arterial_mesenchyme_connectivity.csv -
    which Mesenchyme cluster the Arterial-Endocardium subtype (13.6) connects
    to, using the velocity model's own kNN connectivity graph
    (obsp['connectivities']) rather than 2D FDG proximity; target cluster
    identified dynamically (highest total connectivity), not hardcoded
  - arterial_connectivity/{group}_arterial_mesenchyme_connectivity_barplot.png
  - arterial_connectivity/{group}_arterial_target_cluster_highlight_{fdg,umap}.png

Author: Melanie Smith
Email:  melanie.smith@adelaide.edu.au
Date:   2026-09-23

Usage:
    python 13.7-mesenchymal_cluster_characterisation.py

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

cluster_column = "res_0.5"
min_cluster_size = 20  # exclude tiny stray clusters from the DE / characterisation

#-----------------------------------------------------------------------------------
# Per-group processing
#-----------------------------------------------------------------------------------

for group_ID in group_names:

    print(f"\n=== Processing: {group_ID} ===")

    try:
        full_h5ad_path = os.path.join(out_dir, group_ID, f"{group_ID}_spliced_endo_mesen_subset.h5ad")
        adata_full = sc.read_h5ad(full_h5ad_path)

        mes_dir = os.path.join(out_dir, group_ID, "mesenchymal_subtypes", "cluster_characterisation")
        os.makedirs(mes_dir, exist_ok=True)

        mes = adata_full[adata_full.obs["new_celltypes"] == "Mesenchyme"].copy()

        cluster_sizes = mes.obs[cluster_column].value_counts()
        real_clusters = cluster_sizes[cluster_sizes >= min_cluster_size].index.tolist()
        excluded = cluster_sizes[cluster_sizes < min_cluster_size]
        print(f"{cluster_column} Mesenchyme clusters: {cluster_sizes.to_dict()}")
        if len(excluded):
            print(f"Excluding tiny stray clusters (< {min_cluster_size} cells): {excluded.to_dict()}")

        mes_real = mes[mes.obs[cluster_column].isin(real_clusters)].copy()
        mes_real.obs[cluster_column] = mes_real.obs[cluster_column].cat.remove_unused_categories()

        # ------------------------------------------------------------------
        # One-vs-rest Wilcoxon DE across all real clusters at once
        # ------------------------------------------------------------------
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
            print(f"\n{group_ID} — cluster {cl} ({cluster_sizes[cl]} cells) — top 10 markers:")
            print(top.head(10)[["names", "logfoldchanges", "pvals_adj", "pct_nz_group"]].to_string(index=False))

        marker_df = pd.concat(marker_rows, ignore_index=True)
        marker_df.to_csv(os.path.join(mes_dir, f"{group_ID}_mesenchyme_cluster_markers.csv"), index=False)

        # ------------------------------------------------------------------
        # Summary dotplot: top 5 genes per cluster, across all clusters
        # ------------------------------------------------------------------
        sc.pl.rank_genes_groups_dotplot(
            mes_real, groupby=cluster_column, n_genes=5, show=False,
            title=f"{group_ID} — Mesenchyme {cluster_column} cluster markers"
        )
        plt.savefig(
            os.path.join(mes_dir, f"{group_ID}_mesenchyme_cluster_marker_dotplot.png"),
            dpi=300, bbox_inches="tight"
        )
        plt.close()
        print("Saved marker table and dotplot")

        # ------------------------------------------------------------------
        # FDG / UMAP coloured by res_0.5 cluster, Mesenchyme only,
        # Endocardium shown in grey for spatial context
        # ------------------------------------------------------------------
        velocity_h5ad_path = os.path.join(
            out_dir, group_ID, "velocity_output_dynamical", paga_group,
            f"{group_ID}_velocity_dynamical_{pcs}_{neighbours}_{count_level}_{paga_group}.h5ad"
        )
        adata = sc.read_h5ad(velocity_h5ad_path)
        adata.obs["mesenchyme_cluster_plot"] = np.where(
            adata.obs["new_celltypes"] == "Mesenchyme",
            adata.obs[cluster_column].astype(str),
            "other"
        )
        # keep only the real clusters as distinct categories; anything else (Endocardium,
        # or a tiny stray Mesenchyme cluster) is grouped as "other"
        adata.obs["mesenchyme_cluster_plot"] = np.where(
            adata.obs["mesenchyme_cluster_plot"].isin(real_clusters),
            adata.obs["mesenchyme_cluster_plot"], "other"
        )
        adata.obs["mesenchyme_cluster_plot"] = pd.Categorical(adata.obs["mesenchyme_cluster_plot"])

        cmap = plt.get_cmap("tab10")
        palette = {}
        for i, cl in enumerate(real_clusters):
            palette[cl] = cmap(i % 10)
        palette["other"] = (0.85, 0.85, 0.85, 1.0)
        colours = [palette[c] for c in adata.obs["mesenchyme_cluster_plot"].cat.categories]

        for basis, basis_label in [("draw_graph_fa", "fdg"), ("umap_spliced", "umap")]:
            sc.pl.embedding(
                adata, basis=basis, color="mesenchyme_cluster_plot", palette=colours,
                title=f"{group_ID} — Mesenchyme {cluster_column} clusters — {basis_label}", show=False
            )
            plt.savefig(
                os.path.join(mes_dir, f"{group_ID}_mesenchyme_clusters_{basis_label}.png"),
                dpi=300, bbox_inches="tight"
            )
            plt.close()

        print("Saved cluster highlight plots")

        # ------------------------------------------------------------------
        # Which Mesenchyme cluster does the Arterial-Endocardium subtype
        # (13.6) actually connect to? Uses the velocity model's own kNN
        # connectivity graph (obsp['connectivities']) rather than eyeballing
        # 2D FDG proximity, which can be a stochastic layout artifact. The
        # target cluster is identified dynamically (highest total connectivity)
        # rather than hardcoded, since cluster IDs are arbitrary per genotype.
        # ------------------------------------------------------------------
        conn_dir = os.path.join(out_dir, group_ID, "mesenchymal_subtypes", "arterial_connectivity")
        os.makedirs(conn_dir, exist_ok=True)

        subtype_csv = os.path.join(out_dir, group_ID, "endocardial_subtypes",
                                    f"{group_ID}_endocardial_subtype_assignments.csv")
        subtype_df = pd.read_csv(subtype_csv, index_col=0)
        adata.obs["endocardial_subtype"] = subtype_df.loc[adata.obs_names, "endocardial_subtype"].values

        arterial_mask = (adata.obs["endocardial_subtype"] == "Arterial").values
        mesenchyme_mask = (adata.obs["new_celltypes"] == "Mesenchyme").values

        conn = adata.obsp["connectivities"]
        arterial_to_mesenchyme = conn[arterial_mask][:, mesenchyme_mask]
        mes_cluster_of_cell = adata.obs.loc[mesenchyme_mask, cluster_column].astype(str).values
        conn_weight_per_mes_cell = np.asarray(arterial_to_mesenchyme.sum(axis=0)).flatten()

        conn_df = pd.DataFrame({"mes_cluster": mes_cluster_of_cell, "conn_weight": conn_weight_per_mes_cell})
        cluster_sizes_all = adata.obs.loc[mesenchyme_mask, cluster_column].astype(str).value_counts()
        conn_summary = conn_df.groupby("mes_cluster")["conn_weight"].agg(["sum", "count"])
        conn_summary["mean_conn_per_cell"] = conn_summary["sum"] / conn_summary["count"]
        conn_summary["cluster_size"] = cluster_sizes_all
        conn_summary = conn_summary.sort_values("sum", ascending=False)
        conn_summary.to_csv(os.path.join(conn_dir, f"{group_ID}_arterial_mesenchyme_connectivity.csv"))

        target_cluster = conn_summary.index[0]
        print(f"Arterial-Endocardium connects most strongly to Mesenchyme cluster {target_cluster} "
              f"(total connectivity {conn_summary['sum'].iloc[0]:.2f}, "
              f"next-highest {conn_summary['sum'].iloc[1]:.2f})")
        print(conn_summary)

        # Bar plot of total connectivity by cluster
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.bar(conn_summary.index.astype(str), conn_summary["sum"], color="#E41A1C", alpha=0.8)
        ax.set_xlabel(f"Mesenchyme {cluster_column} cluster")
        ax.set_ylabel("Total connectivity from Arterial-Endocardium")
        ax.set_title(f"{group_ID} — Arterial-Endocardium connectivity by Mesenchyme cluster")
        fig.tight_layout()
        fig.savefig(
            os.path.join(conn_dir, f"{group_ID}_arterial_mesenchyme_connectivity_barplot.png"),
            dpi=300, bbox_inches="tight"
        )
        plt.close(fig)

        # FDG/UMAP highlight: Arterial-Endocardium cells + target Mesenchyme cluster
        adata.obs["arterial_target_plot"] = "other"
        adata.obs.loc[mesenchyme_mask & (adata.obs[cluster_column].astype(str) == target_cluster),
                      "arterial_target_plot"] = f"Mesenchyme cluster {target_cluster} (target)"
        adata.obs.loc[arterial_mask, "arterial_target_plot"] = "Arterial-Endocardium"
        adata.obs["arterial_target_plot"] = pd.Categorical(adata.obs["arterial_target_plot"])

        highlight_palette = {
            "other": (0.85, 0.85, 0.85, 1.0),
            f"Mesenchyme cluster {target_cluster} (target)": "#377EB8",
            "Arterial-Endocardium": "#E41A1C",
        }
        colours = [highlight_palette[c] for c in adata.obs["arterial_target_plot"].cat.categories]

        for basis, basis_label in [("draw_graph_fa", "fdg"), ("umap_spliced", "umap")]:
            sc.pl.embedding(
                adata, basis=basis, color="arterial_target_plot", palette=colours,
                title=f"{group_ID} — Arterial-Endocardium & its target Mesenchyme cluster — {basis_label}",
                show=False
            )
            plt.savefig(
                os.path.join(conn_dir, f"{group_ID}_arterial_target_cluster_highlight_{basis_label}.png"),
                dpi=300, bbox_inches="tight"
            )
            plt.close()

        print("Saved Arterial-Mesenchyme connectivity outputs")
        print(f"Finished: {group_ID}")

    except Exception as e:
        plt.close("all")
        warnings.warn(f"{group_ID} — mesenchymal cluster characterisation failed: {e}")
