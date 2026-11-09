#!/usr/bin/env python3
"""
Additional plots requested by Sophie after the 14.x/15.x review meeting.
Not a new analysis - reads the already-computed 14.4 final velocity object
and makes extra views of it, with no re-running of the dynamical model fit.

This round (task 1 of 5 from the meeting notes, plus a confidence check):
  - FDG velocity embedding, split into separate panels by genotype
    (control vs Mib1_KO), to directly compare velocity structure between
    genotypes
  - velocity confidence (scv.tl.velocity_confidence), on both bases, as a
    scale-independent complement to eyeballing arrow length in the velocity
    plots - arrow length in scv.pl.velocity_embedding is a real, per-cell
    data-driven vector, but it's uniformly rescaled by a fixed display
    constant (arrow_length) for legibility, so it's only directly comparable
    within one plot. velocity_confidence gives an explicit per-cell number
    instead, independent of any plotting scale.
  - task 2: UMAP feature plots for NCC marker genes Plxna2 and Sema3c, split
    by sample and by timepoint (stage). Uses sc.pl.embedding (not
    scv.pl.scatter, which has no vmin/vmax) with an explicit, shared
    vmin/vmax per gene computed once across all cells - otherwise each
    per-panel call would independently rescale its own colour range to its
    own subset's min/max, the same bug fixed for the EGFP split-by-sample
    plot in 14.1/15.1.
  - task 3: latent time (UMAP + FDG), restricted to Mesenchyme cells only -
    a single subset per basis, so (unlike tasks 1/2) no panel-split colour
    scale issue to guard against.

Title/legend fonts are set larger throughout this script so the plots stay
legible once pasted into the lab book at a smaller size.

Author: Melanie Smith
Email:  melanie.smith@adelaide.edu.au
Date:   2026-10-15

Usage:
    python 14.6-additional_velocity_plots_for_sophie.py
"""

#-----------------------------------------------------------------------------------
#   Import core libraries
#-----------------------------------------------------------------------------------

import scvelo as scv
import scanpy as sc
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

title_fontsize = 20
legend_fontsize = 14

project_dir = "/home/melanie-smith/workDir/sophieWiszniak/20251104_sophieWiszniak_ncc_pa"
out_dir     = os.path.join(project_dir, "outDir", "14-rna_velocity")
velocity_out_dir = os.path.join(out_dir, "velocity_output_dynamical", paga_group)

h5ad_path = os.path.join(
    velocity_out_dir,
    f"all_samples_velocity_dynamical_{pcs}_{neighbours}_{count_level}_{paga_group}.h5ad"
)

print("\n=== Processing: all_samples ===")
adata = sc.read_h5ad(h5ad_path)

#-----------------------------------------------------------------------------------
# Task 1 (Sophie, meeting notes): velocity embedding split by genotype and by
# sample, on both the FDG and UMAP bases. Velocity graph/embedding was
# already computed on the full merged object in 14.4 - subsetting here only
# selects which already-computed cells/vectors to draw, it does not
# recompute the velocity graph per panel.
#-----------------------------------------------------------------------------------

def plot_velocity_split(adata, basis, basis_title, split_col, split_values, nrows, ncols, figsize, filename):
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    axes_flat = axes.flat if hasattr(axes, "flat") else [axes]
    for i, (ax, value) in enumerate(zip(axes_flat, split_values)):
        scv.pl.velocity_embedding(
            adata[adata.obs[split_col] == value],
            basis=basis,
            color="new_celltypes",
            arrow_length=3,
            arrow_size=2,
            title=value,
            fontsize=title_fontsize,
            legend_fontsize=legend_fontsize,
            legend_loc="right margin" if i == len(split_values) - 1 else "none",
            ax=ax,
            show=False,
        )
    fig.suptitle(
        f"{basis_title} velocity embedding by {split_col} - {count_level} - all_samples",
        fontsize=title_fontsize
    )
    fig.savefig(os.path.join(velocity_out_dir, filename), dpi=300, bbox_inches="tight")
    plt.close(fig)


genotypes = ["control", "Mib1_KO"]
samples = ["e11_control", "e11_ko", "e12_control", "e12_ko"]

plot_velocity_split(
    adata, "draw_graph_fa", "FDG", "genotype", genotypes, 1, 2, (16, 7),
    f"all_samples_velocity_fdg_{pcs}_{neighbours}_{count_level}_{paga_group}_by_genotype.png"
)
plot_velocity_split(
    adata, "draw_graph_fa", "FDG", "orig.ident", samples, 2, 2, (16, 14),
    f"all_samples_velocity_fdg_{pcs}_{neighbours}_{count_level}_{paga_group}_by_sample.png"
)
plot_velocity_split(
    adata, "umap_spliced", "UMAP", "genotype", genotypes, 1, 2, (16, 7),
    f"all_samples_velocity_umap_{count_level}_{paga_group}_by_genotype.png"
)
plot_velocity_split(
    adata, "umap_spliced", "UMAP", "orig.ident", samples, 2, 2, (16, 14),
    f"all_samples_velocity_umap_{count_level}_{paga_group}_by_sample.png"
)

#-----------------------------------------------------------------------------------
# Velocity confidence - scale-independent complement to arrow length
#-----------------------------------------------------------------------------------

scv.tl.velocity_confidence(adata)

for basis, basis_label, extra in [
    ("umap_spliced", "umap", ""),
    ("draw_graph_fa", "fdg", f"_{pcs}_{neighbours}"),
]:
    scv.pl.scatter(
        adata,
        color="velocity_confidence",
        basis=basis,
        cmap="coolwarm",
        title=f"Velocity confidence on {basis_label.upper()} - all_samples",
        fontsize=title_fontsize,
        legend_fontsize=legend_fontsize,
        show=False
    )
    plt.savefig(
        os.path.join(
            velocity_out_dir,
            f"all_samples_qc_velocity_confidence_{basis_label}{extra}_{count_level}_{paga_group}.png"
        ),
        dpi=300, bbox_inches="tight"
    )
    plt.close()

#-----------------------------------------------------------------------------------
# Task 2 (Sophie, meeting notes): UMAP feature plots for Plxna2/Sema3c,
# split by sample and by timepoint. Shared vmin/vmax computed once per gene
# across all cells, then passed explicitly to every panel - sc.pl.embedding
# (unlike scv.pl.scatter) accepts vmin/vmax directly.
#-----------------------------------------------------------------------------------

stages = ["e11", "e12"]
ncc_genes = ["Plxna2", "Sema3c"]

def plot_expression_split(adata, basis, basis_title, gene, vmin, vmax, split_col, split_values, nrows, ncols, figsize, filename):
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    axes_flat = axes.flat if hasattr(axes, "flat") else [axes]
    for ax, value in zip(axes_flat, split_values):
        sc.pl.embedding(
            adata[adata.obs[split_col] == value],
            basis=basis,
            color=gene,
            vmin=vmin,
            vmax=vmax,
            title=value,
            size=18,
            frameon=False,
            ax=ax,
            show=False,
        )
        ax.set_title(value, fontsize=title_fontsize)
    fig.suptitle(f"{gene} expression ({basis_title}) by {split_col} - all_samples", fontsize=title_fontsize)
    fig.savefig(os.path.join(velocity_out_dir, filename), dpi=300, bbox_inches="tight")
    plt.close(fig)

for gene in ncc_genes:
    vmin, vmax = adata[:, gene].X.min(), adata[:, gene].X.max()
    plot_expression_split(
        adata, "umap_spliced", "UMAP", gene, vmin, vmax, "orig.ident", samples, 2, 2, (16, 14),
        f"all_samples_{gene}_umap_{count_level}_{paga_group}_by_sample.png"
    )
    plot_expression_split(
        adata, "umap_spliced", "UMAP", gene, vmin, vmax, "stage", stages, 1, 2, (16, 7),
        f"all_samples_{gene}_umap_{count_level}_{paga_group}_by_timepoint.png"
    )

#-----------------------------------------------------------------------------------
# Task 3 (Sophie, meeting notes): latent time (UMAP + FDG), Mesenchyme only.
#-----------------------------------------------------------------------------------

mesenchyme_adata = adata[adata.obs["new_celltypes"] == "Mesenchyme"]

for basis, basis_label, extra in [
    ("umap_spliced", "umap", ""),
    ("draw_graph_fa", "fdg", f"_{pcs}_{neighbours}"),
]:
    scv.pl.scatter(
        mesenchyme_adata,
        color="latent_time",
        basis=basis,
        cmap="gnuplot",
        title=f"Latent time on {basis_label.upper()} - Mesenchyme only - all_samples",
        fontsize=title_fontsize,
        legend_fontsize=legend_fontsize,
        show=False
    )
    plt.savefig(
        os.path.join(
            velocity_out_dir,
            f"all_samples_mesenchyme_only_qc_latenttime_{basis_label}{extra}_{count_level}_{paga_group}.png"
        ),
        dpi=300, bbox_inches="tight"
    )
    plt.close()

print("Finished: all_samples")
