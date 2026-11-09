#!/usr/bin/env python3
"""
Follow-up plot requested at the 2026-10-26 meeting with Sophie, not added to
the additional_plots_meeting_20261013 lab book (saved locally only, to be
added by hand) - latent time (Endocardium + Mesenchyme object, both cell
types) split into 4 panels by orig.ident, on both the FDG and UMAP bases,
both lineages (14.x cc-retained, 15.x no cc).

Reads the already-fit 14.4/15.4 final velocity .h5ad objects - no
re-running the dynamical model fit, just re-plotting the already-computed
latent_time. Shared vmin/vmax=(0,1) across all 4 panels (latent_time is
already on a 0-1 scale) so panels are directly comparable, same colour-scale
discipline as the by-sample/by-timepoint gene expression plots in 14.6/15.5.

Author: Melanie Smith
Email:  melanie.smith@adelaide.edu.au
Date:   2026-11-02

Usage:
    python 14.15-latenttime_by_sample.py
"""

import scvelo as scv
import scanpy as sc
import matplotlib.pyplot as plt
import os

paga_group = "celltypes"
count_level = "spliced_unspliced"
pcs_value = 30
neighbours_value = 25
pcs = f"pcs{pcs_value}"
neighbours = f"neighbours{neighbours_value}"

title_fontsize = 20
legend_fontsize = 14

project_dir = "/home/melanie-smith/workDir/sophieWiszniak/20251104_sophieWiszniak_ncc_pa"

samples = ["e11_control", "e11_ko", "e12_control", "e12_ko"]

lineages = [
    ("14-rna_velocity", "cc-retained"),
    ("15-rna_velocity", "no cc"),
]

for out_dir_name, lineage_label in lineages:
    print(f"\n=== Processing: {out_dir_name} ({lineage_label}) ===")

    out_dir = os.path.join(project_dir, "outDir", out_dir_name)
    velocity_out_dir = os.path.join(out_dir, "velocity_output_dynamical", paga_group)

    h5ad_path = os.path.join(
        velocity_out_dir,
        f"all_samples_velocity_dynamical_{pcs}_{neighbours}_{count_level}_{paga_group}.h5ad"
    )
    adata = sc.read_h5ad(h5ad_path)

    for basis, basis_label, extra in [
        ("umap_spliced", "umap", ""),
        ("draw_graph_fa", "fdg", f"_{pcs}_{neighbours}"),
    ]:
        fig, axes = plt.subplots(2, 2, figsize=(16, 14))
        for ax, sample in zip(axes.flat, samples):
            scv.pl.scatter(
                adata[adata.obs["orig.ident"] == sample],
                color="latent_time",
                basis=basis,
                cmap="gnuplot",
                vmin=0,
                vmax=1,
                title=sample,
                fontsize=title_fontsize,
                legend_fontsize=legend_fontsize,
                ax=ax,
                show=False,
            )
        fig.suptitle(
            f"Latent time on {basis_label.upper()} by sample - Endocardium + Mesenchyme ({lineage_label}) - all_samples",
            fontsize=title_fontsize
        )
        fig.savefig(
            os.path.join(
                velocity_out_dir,
                f"all_samples_qc_latenttime_{basis_label}{extra}_{count_level}_{paga_group}_by_sample.png"
            ),
            dpi=300, bbox_inches="tight"
        )
        plt.close(fig)

    print(f"Finished: {out_dir_name}")
