#!/usr/bin/env python3
"""
SCVELO
=================

Runs RNA velocity in dynamical mode on the PAGA/FDG-processed spliced/unspliced
Endocardium + Mesenchyme object (from 14.3-run_scanpy_paga_fdg_cell.py), merging
in the raw Velocyto .loom spliced/unspliced/ambiguous counts.

Mirrors 12.4/13.4-run_rnaVelocity_dynamical_fdg_celltypes.py, updated to:
  - point at the 14.3-generated PAGA/FDG .h5ad object
    (outDir/14-rna_velocity/scanpy_output/celltypes/...) - the all-four-
    samples-merged, EGFP-lineage-tracing counterpart of 12.3/13.3's
    genotype-pooled objects
  - no outer loop over groups - a single fixed (pcs=30, neighbours=25)
    combination is used directly on the one merged object
  - merge in .loom spliced/unspliced/ambiguous counts per constituent sample
    (identified via 'orig.ident'), same as 12.4/13.4, just four samples
    merged here instead of two per group
  - use 'seurat_clusters' for the cluster QC plot (this lineage, like 12.x,
    retains cell-cycle genes, so 'seurat_clusters' is the direct analog
    carried through from the source object - not the 'seurat_clusters_no_cc'
    used in the 13.x no-cc lineage)
  - add EGFP-specific plots: a velocity embedding (arrows) coloured by
    EGFP_UMI_log1p on both FDG and UMAP, using the same grey90/gold/red
    colour scale established in 14.1/14.2 (Spectral-style palettes have no
    neutral low end, which misleadingly colours zero-value cells for a
    sparse feature like this), plus simple QC scatter plots for
    EGFP_positive on both bases

Author: Melanie Smith
Email:  melanie.smith@adelaide.edu.au
Date:   2026-10-07

Usage:
    python 14.4-run_rnaVelocity_dynamical_fdg_celltypes.py

Notes:
    - Main libraries used: scvelo, scanpy, anndata, matplotlib
    - Python version: >= 3.9 recommended
"""

#-----------------------------------------------------------------------------------
#   Import core libraries
#-----------------------------------------------------------------------------------

import scvelo as scv
import scanpy as sc # single-cell analysis toolkit
import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import os
import warnings

#-----------------------------------------------------------------------------------
#   Params
#-----------------------------------------------------------------------------------

paga_group = "celltypes"
count_level = "spliced_unspliced"

# Fixed (pcs, neighbours) combination, chosen after inspecting the 14.3 FDG
# sweep plots, consistent with the (pcs=30, neighbours=25) used for 12.4/13.4
pcs_value = 30
neighbours_value = 25
pcs = f"pcs{pcs_value}"
neighbours = f"neighbours{neighbours_value}"

project_dir = "/home/melanie-smith/workDir/sophieWiszniak/20251104_sophieWiszniak_ncc_pa"
out_dir     = os.path.join(project_dir, "outDir", "14-rna_velocity")

sample_names = ["e11_control", "e11_ko", "e12_control", "e12_ko"]

# Named vector mapping each original sample to its Velocyto loom column prefix
# (same mapping used throughout the R/Python pipeline for this project)
loom_colnames = {
    "e11_control": "possorted_genome_bam_5Q1UF:",
    "e11_ko":      "possorted_genome_bam_LZRTJ:",
    "e12_control": "possorted_genome_bam_YIL7N:",
    "e12_ko":      "possorted_genome_bam_0YX73:",
}

# grey90 -> gold -> red: matches the colour scale established in 14.1/14.2 for
# EGFP (a Spectral-style palette has no neutral low end, so zero-value cells
# for this sparse feature still render as a solid, misleading colour)
egfp_cmap = LinearSegmentedColormap.from_list("egfp", ["#e5e5e5", "gold", "red"])

#-----------------------------------------------------------------------------------
# Per-sample loom merge helper
#-----------------------------------------------------------------------------------

def merge_sample_loom(adata_all, sample_ID, project_dir):
    """
    Subset adata_all to sample_ID's cells (via 'orig.ident'), match those
    cells against sample_ID's own raw .loom file on the bare barcode (safe -
    no collisions within a single sample), merge in the spliced/unspliced/
    ambiguous layers, then restore each cell's original, dataset-unique
    barcode name before handing back to the caller for concatenation across
    samples.
    """

    adata_sample = adata_all[adata_all.obs["orig.ident"] == sample_ID].copy()
    original_names = adata_sample.obs_names.copy()

    # Bare barcode for matching against this sample's own .loom file -
    # splitting on "_" and taking the first segment keeps the "-1" suffix
    # (which comes before any Seurat-added multi-sample "_N" suffix)
    bare_names = adata_sample.obs_names.str.split("_").str[0]
    adata_sample.obs_names = bare_names

    loom_file = os.path.join(
        project_dir, "outDir", "09-rna_velocity", sample_ID, "velocyto_output",
        f"{loom_colnames[sample_ID].rstrip(':')}.loom"
    )
    vlm = sc.read_loom(loom_file)

    # Remove sample prefix, trailing x, add -1 suffix
    vlm.obs_names = vlm.obs_names.str.split(":").str[1]
    vlm.obs_names = vlm.obs_names.str.rstrip("x")
    vlm.obs_names = vlm.obs_names + "-1"

    common = adata_sample.obs_names.intersection(vlm.obs_names)
    pct_matched = 100 * len(common) / len(adata_sample.obs_names)
    print(f"  {sample_ID}: adata={len(adata_sample)}, loom={len(vlm)}, common={len(common)} ({pct_matched:.1f}%)")

    if pct_matched < 95:
        warnings.warn(
            f"{sample_ID} — only {pct_matched:.1f}% of cells matched to the .loom file"
        )

    adata_sample = adata_sample[common].copy()
    vlm = vlm[common].copy()
    vlm = vlm[adata_sample.obs_names].copy()

    merged_sample = scv.utils.merge(adata_sample, vlm)

    # Restore original, dataset-unique barcode names (mapped through the same
    # bare-barcode intersection/order used above) so the four samples' merged
    # subsets can be concatenated without name collisions.
    name_map = pd.Series(original_names.values, index=bare_names.values)
    merged_sample.obs_names = name_map.loc[merged_sample.obs_names].values

    return merged_sample

#-----------------------------------------------------------------------------------
# Run RNA velocity (dynamical mode)
#-----------------------------------------------------------------------------------

print("\n=== Processing: all_samples ===")

try:
    # set the .h5ad data input file address - the 14.3-generated PAGA/FDG object
    h5ad_file_input = os.path.join(
        out_dir, "scanpy_output", paga_group,
        f"all_samples_fdg_{pcs}_{neighbours}_counts_{count_level}_{paga_group}.h5ad"
    )

    # set the output directory
    velocity_out_dir = os.path.join(out_dir, "velocity_output_dynamical", paga_group)
    os.makedirs(velocity_out_dir, exist_ok=True)

    #-----------------------------------------------------------------------------------
    # Import the scv matrix and merge in each constituent sample's .loom data
    #-----------------------------------------------------------------------------------

    adata_all = sc.read_h5ad(h5ad_file_input)
    print("h5ad example:", adata_all.obs_names[0])

    print("Merging .loom data per constituent sample:")
    merged_parts = [
        merge_sample_loom(adata_all, sample_ID, project_dir)
        for sample_ID in sample_names
    ]

    adata = ad.concat(merged_parts)

    # safety check after merging
    print(adata.layers.keys())
    print("merged total:", len(adata))

    #-----------------------------------------------------------------------------------
    # let's take a look at the proportion of spliced and unspliced reads
    #-----------------------------------------------------------------------------------

    scv.pl.proportions(adata, show=False)
    plt.savefig(
        os.path.join(velocity_out_dir, f"all_samples_proportions_{pcs}_{neighbours}_{count_level}_{paga_group}.png"),
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()

    #-----------------------------------------------------------------------------------
    # now we need to normalise the spliced data
    #-----------------------------------------------------------------------------------

    # ensure previous moments are removed (I don't plan to re-run but just in case)
    adata.layers.pop("Ms", None)
    adata.layers.pop("Mu", None)
    adata.uns.pop("neighbors", None)

    # filter and normalise
    scv.pp.filter_and_normalize(adata)

    #-----------------------------------------------------------------------------------
    # calculate the moments (i.e. means and uncentered variances) computed among nearest neighbors in PCA space
    #-----------------------------------------------------------------------------------
    scv.pp.moments(adata)

    #-----------------------------------------------------------------------------------
    # compute the dynamical model and velocities
    #-----------------------------------------------------------------------------------

    # find the top genes, quicker and less noisy than using all genes
    scv.pp.filter_genes(adata, min_shared_counts=20)
    scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)

    scv.tl.recover_dynamics(adata, n_jobs=4, show_progress_bar=False)          # this can take time
    scv.tl.velocity(adata, mode='dynamical')

    #-----------------------------------------------------------------------------------
    # compute the velocity graph
    #-----------------------------------------------------------------------------------

    scv.tl.velocity_graph(adata)

    # Anchor latent_time's direction to a single, fixed root cell rather
    # than letting the model pick among Endocardium candidates via its own
    # (ambiguous) fitted gene-time sum. Choose the Endocardium cell most
    # deeply embedded within the Endocardium population - highest
    # connectivity-graph density of other Endocardium cells - as a robust,
    # reproducible representative root, then pass it as a literal index via
    # adata.uns['iroot'].
    root_mask = (adata.obs['new_celltypes'] == 'Endocardium').astype(float).values
    root_density = adata.obsp['connectivities'].dot(root_mask)
    root_density[root_mask == 0] = -np.inf  # restrict the choice to Endocardium cells themselves
    adata.uns['iroot'] = int(np.argmax(root_density))

    scv.tl.latent_time(adata, root_key='iroot')

    #-----------------------------------------------------------------------------------
    # plot velocity arrows on the UMAP structure (umap_spliced - no gene-level
    # 'umap' exists in this object)
    #-----------------------------------------------------------------------------------

    scv.pl.velocity_embedding(
        adata,
        title=f"umap_spliced projection - dynamical - {count_level} - all_samples",
        basis='umap_spliced',
        color='new_celltypes',
        arrow_length=3,
        arrow_size=2,
        dpi=120,
        show=False,
        legend_loc="bottom right"
    )
    plt.savefig(os.path.join(velocity_out_dir, f"all_samples_velocity_umap_{count_level}_{paga_group}.png"),
                dpi=300,
                bbox_inches="tight")
    plt.close()

    # QC plots
    scv.pl.scatter(adata, color='new_celltypes', basis='umap_spliced', title='Celltypes on UMAP', show=False)
    plt.savefig(os.path.join(velocity_out_dir, f"all_samples_qc_celltypes_umap_{count_level}_{paga_group}.png"),
                dpi=300,
                bbox_inches="tight")
    plt.close()

    scv.pl.scatter(adata, color='seurat_clusters', basis='umap_spliced', title='Clusters on UMAP', show=False)
    plt.savefig(os.path.join(velocity_out_dir, f"all_samples_qc_clusters_umap_{count_level}_{paga_group}.png"),
                dpi=300,
                bbox_inches="tight")
    plt.close()

    scv.pl.scatter(adata, color='latent_time', basis='umap_spliced', cmap='gnuplot', title='Latent time on UMAP', show=False)
    plt.savefig(os.path.join(velocity_out_dir, f"all_samples_qc_latenttime_umap_{count_level}_{paga_group}.png"),
                dpi=300,
                bbox_inches="tight")
    plt.close()

    #-----------------------------------------------------------------------------------
    # plot velocity arrows on the FDG structure
    #-----------------------------------------------------------------------------------

    scv.pl.velocity_embedding(
        adata,
        title=f"PCs= {pcs}, Neighbours= {neighbours} - ({paga_group}) - {count_level} - dynamical - all_samples",
        basis='draw_graph_fa',
        color='new_celltypes',
        arrow_length=3,
        arrow_size=2,
        dpi=120,
        show=False,
        legend_loc="bottom right"
    )
    plt.savefig(os.path.join(velocity_out_dir, f"all_samples_velocity_fdg_{pcs}_{neighbours}_{count_level}_{paga_group}.png"),
                dpi=300,
                bbox_inches="tight")
    plt.close()

    # QC plots
    scv.pl.scatter(adata, color='new_celltypes', basis='draw_graph_fa', title='Celltypes on FDG', show=False)
    plt.savefig(os.path.join(velocity_out_dir, f"all_samples_qc_celltypes_fdg_{pcs}_{neighbours}_{count_level}_{paga_group}.png"),
                dpi=300,
                bbox_inches="tight")
    plt.close()

    scv.pl.scatter(adata, color='seurat_clusters', basis='draw_graph_fa', title='Clusters on FDG', show=False)
    plt.savefig(os.path.join(velocity_out_dir, f"all_samples_qc_clusters_fdg_{pcs}_{neighbours}_{count_level}_{paga_group}.png"),
                dpi=300,
                bbox_inches="tight")
    plt.close()

    scv.pl.scatter(adata, color='latent_time', basis='draw_graph_fa', cmap='gnuplot', title='Latent time on FDG', show=False)
    plt.savefig(os.path.join(velocity_out_dir, f"all_samples_qc_latenttime_fdg_{pcs}_{neighbours}_{count_level}_{paga_group}.png"),
                dpi=300,
                bbox_inches="tight")
    plt.close()

    #-----------------------------------------------------------------------------------
    # EGFP - velocity embedding (arrows) coloured by EGFP_UMI, plus QC scatter
    # for EGFP_positive, on both bases
    #-----------------------------------------------------------------------------------

    for basis, basis_label in [("umap_spliced", "umap"), ("draw_graph_fa", "fdg")]:

        scv.pl.velocity_embedding(
            adata,
            title=f"EGFP_UMI - {basis} - dynamical - {count_level} - all_samples",
            basis=basis,
            color='EGFP_UMI_log1p',
            cmap=egfp_cmap,
            arrow_length=3,
            arrow_size=2,
            dpi=120,
            show=False
        )
        plt.savefig(
            os.path.join(velocity_out_dir, f"all_samples_velocity_{basis_label}_{count_level}_{paga_group}_EGFP_UMI.png"),
            dpi=300, bbox_inches="tight"
        )
        plt.close()

        scv.pl.scatter(adata, color='EGFP_positive', basis=basis,
                        title=f"EGFP positive on {basis_label.upper()}", show=False)
        plt.savefig(
            os.path.join(velocity_out_dir, f"all_samples_qc_egfp_positive_{basis_label}_{count_level}_{paga_group}.png"),
            dpi=300, bbox_inches="tight"
        )
        plt.close()

    #-----------------------------------------------------------------------------------
    # write out the velocity-processed object for downstream use
    #-----------------------------------------------------------------------------------

    adata.write_h5ad(
        os.path.join(velocity_out_dir, f"all_samples_velocity_dynamical_{pcs}_{neighbours}_{count_level}_{paga_group}.h5ad"),
        compression="gzip"
    )

    print("Finished: all_samples")

except Exception as e:
    plt.close("all")  # make sure no partially-drawn figure is left open
    warnings.warn(f"all_samples — velocity run failed: {e}")
