#!/usr/bin/env python3
"""
SCVELO
=================

Runs RNA velocity in dynamical mode on the PAGA/FDG-processed spliced/unspliced
Endocardium + Mesenchyme genotype-group objects (from
12.3-run_scanpy_paga_fdg_cell.py), merging in the raw Velocyto .loom
spliced/unspliced/ambiguous counts.

Mirrors 11.4-run_rnaVelocity_dynamical_fdg_celltypes.py, updated to:
  - point at the 12.3-generated PAGA/FDG .h5ad objects
    (outDir/12-rna_velocity/{group_ID}/scanpy_output/celltypes/...)
  - loop over the two genotype groups ("control", "ko") instead of the four
    individual samples, using a single fixed (pcs=30, neighbours=25)
    combination for both groups
  - merge in .loom spliced/unspliced/ambiguous counts per constituent sample,
    not per group: each group's .h5ad pools two original samples
    (e.g. control = e11_control + e12_control), and raw 10x barcodes collide
    between samples once stripped to their bare form (confirmed empirically -
    3 collisions among 7599 control cells). Barcode matching against each
    sample's own .loom file is therefore done on a per-sample subset
    (identified via 'orig.ident', not by inferring from barcode-suffix
    patterns), using each cell's *original* group-level barcode as the
    identifier once merged, before the two merged sample subsets are
    concatenated back into one group-level object.
  - use 'seurat_clusters' for the cluster QC plot instead of
    'full_dataset_clusters_no_cc' - the 12.x lineage never removed cell-cycle
    genes, so there is no "_no_cc" variant; 'seurat_clusters' is the direct
    analog already carried through from the source object

Author: Melanie Smith
Email:  melanie.smith@adelaide.edu.au
Date:   2026-09-22

Usage:
    python 12.4-run_rnaVelocity_dynamical_fdg_celltypes.py

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
import numpy as np
import os
import warnings

#-----------------------------------------------------------------------------------
#   Params
#-----------------------------------------------------------------------------------

paga_group = "celltypes"
count_level = "spliced_unspliced"

group_names = ["control", "ko"]

# Fixed (pcs, neighbours) combination used for both groups, chosen as a single
# representative combination after inspecting the 12.3 FDG sweep plots, rather
# than repeating the expensive dynamical model fit across all 16 combinations.
pcs_value = 30
neighbours_value = 25
pcs = f"pcs{pcs_value}"
neighbours = f"neighbours{neighbours_value}"

project_dir = "/home/melanie-smith/workDir/sophieWiszniak/20251104_sophieWiszniak_ncc_pa"

# Each genotype group pools two original samples (time points E11 + E12)
group_sample_map = {
    "control": ["e11_control", "e12_control"],
    "ko":      ["e11_ko", "e12_ko"],
}

# Named vector mapping each original sample to its Velocyto loom column prefix
# (same mapping used throughout the R/Python pipeline for this project)
loom_colnames = {
    "e11_control": "possorted_genome_bam_5Q1UF:",
    "e11_ko":      "possorted_genome_bam_LZRTJ:",
    "e12_control": "possorted_genome_bam_YIL7N:",
    "e12_ko":      "possorted_genome_bam_0YX73:",
}

#-----------------------------------------------------------------------------------
# Per-sample loom merge helper
#-----------------------------------------------------------------------------------

def merge_sample_loom(adata_group, sample_ID, project_dir):
    """
    Subset adata_group to sample_ID's cells (via 'orig.ident'), match those
    cells against sample_ID's own raw .loom file on the bare barcode (safe -
    no collisions within a single sample), merge in the spliced/unspliced/
    ambiguous layers, then restore each cell's original, group-unique barcode
    name before handing back to the caller for concatenation across samples.
    """

    adata_sample = adata_group[adata_group.obs["orig.ident"] == sample_ID].copy()
    original_names = adata_sample.obs_names.copy()

    # Bare barcode for matching against this sample's own .loom file
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

    # Restore original, group-unique barcode names (mapped through the same
    # bare-barcode intersection/order used above) so the two samples' merged
    # subsets can be concatenated without name collisions.
    name_map = pd.Series(original_names.values, index=bare_names.values)
    merged_sample.obs_names = name_map.loc[merged_sample.obs_names].values

    return merged_sample

#-----------------------------------------------------------------------------------
# Run RNA velocity (dynamical mode) for each genotype group
#-----------------------------------------------------------------------------------

for group_ID in group_names:

    print(f"\n=== Processing: {group_ID} ===")

    try:
        # set the .h5ad data input file address - the 12.3-generated PAGA/FDG object
        h5ad_file_input = os.path.join(
            project_dir, "outDir", "12-rna_velocity", group_ID, "scanpy_output", paga_group,
            f"{group_ID}_fdg_{pcs}_{neighbours}_counts_{count_level}_{paga_group}.h5ad"
        )

        # set the output directory
        out_dir = os.path.join(
            project_dir, "outDir", "12-rna_velocity", group_ID, "velocity_output_dynamical", paga_group
        )
        os.makedirs(out_dir, exist_ok=True)

        #-----------------------------------------------------------------------------------
        # Import the scv matrix and merge in each constituent sample's .loom data
        #-----------------------------------------------------------------------------------

        adata_group = sc.read_h5ad(h5ad_file_input)
        print("h5ad example:", adata_group.obs_names[0])

        print(f"Merging .loom data per constituent sample for {group_ID}:")
        merged_parts = [
            merge_sample_loom(adata_group, sample_ID, project_dir)
            for sample_ID in group_sample_map[group_ID]
        ]

        adata = ad.concat(merged_parts)

        # safety check after merging
        print(adata.layers.keys())
        print("merged group total:", len(adata))

        #-----------------------------------------------------------------------------------
        # let's take a look at the proportion of spliced and unspliced reads
        #-----------------------------------------------------------------------------------

        scv.pl.proportions(adata, show=False)
        plt.savefig(
            os.path.join(out_dir, f"{group_ID}_proportions_{pcs}_{neighbours}_{count_level}_{paga_group}.png"),
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
        # (ambiguous) fitted gene-time sum, which only partially corrected the
        # direction. Choose the Endocardium cell most deeply embedded within
        # the Endocardium population - highest connectivity-graph density of
        # other Endocardium cells - as a robust, reproducible representative
        # root, then pass it as a literal index via adata.uns['iroot'].
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
            title=f"umap_spliced projection - dynamical - {count_level} - {group_ID}",
            basis='umap_spliced',
            color='new_celltypes',
            arrow_length=3,
            arrow_size=2,
            dpi=120,
            show=False,
            legend_loc="bottom right"
        )
        plt.savefig(os.path.join(out_dir, f"{group_ID}_velocity_umap_{count_level}_{paga_group}.png"),
                    dpi=300,
                    bbox_inches="tight")
        plt.close()

        # QC plots
        scv.pl.scatter(adata, color='new_celltypes', basis='umap_spliced', title='Celltypes on UMAP', show=False)
        plt.savefig(os.path.join(out_dir, f"{group_ID}_qc_celltypes_umap_{count_level}_{paga_group}.png"),
                    dpi=300,
                    bbox_inches="tight")
        plt.close()

        scv.pl.scatter(adata, color='seurat_clusters', basis='umap_spliced', title='Clusters on UMAP', show=False)
        plt.savefig(os.path.join(out_dir, f"{group_ID}_qc_clusters_umap_{count_level}_{paga_group}.png"),
                    dpi=300,
                    bbox_inches="tight")
        plt.close()

        scv.pl.scatter(adata, color='latent_time', basis='umap_spliced', cmap='gnuplot', title='Latent time on UMAP', show=False)
        plt.savefig(os.path.join(out_dir, f"{group_ID}_qc_latenttime_umap_{count_level}_{paga_group}.png"),
                    dpi=300,
                    bbox_inches="tight")
        plt.close()

        #-----------------------------------------------------------------------------------
        # plot velocity arrows on the FDG structure
        #-----------------------------------------------------------------------------------

        scv.pl.velocity_embedding(
            adata,
            title=f"PCs= {pcs}, Neighbours= {neighbours} - ({paga_group}) - {count_level} - dynamical - {group_ID}",
            basis='draw_graph_fa',
            color='new_celltypes',
            arrow_length=3,
            arrow_size=2,
            dpi=120,
            show=False,
            legend_loc="bottom right"
        )
        plt.savefig(os.path.join(out_dir, f"{group_ID}_velocity_fdg_{pcs}_{neighbours}_{count_level}_{paga_group}.png"),
                    dpi=300,
                    bbox_inches="tight")
        plt.close()

        # QC plots
        scv.pl.scatter(adata, color='new_celltypes', basis='draw_graph_fa', title='Celltypes on FDG', show=False)
        plt.savefig(os.path.join(out_dir, f"{group_ID}_qc_celltypes_fdg_{pcs}_{neighbours}_{count_level}_{paga_group}.png"),
                    dpi=300,
                    bbox_inches="tight")
        plt.close()

        scv.pl.scatter(adata, color='seurat_clusters', basis='draw_graph_fa', title='Clusters on FDG', show=False)
        plt.savefig(os.path.join(out_dir, f"{group_ID}_qc_clusters_fdg_{pcs}_{neighbours}_{count_level}_{paga_group}.png"),
                    dpi=300,
                    bbox_inches="tight")
        plt.close()

        scv.pl.scatter(adata, color='latent_time', basis='draw_graph_fa', cmap='gnuplot', title='Latent time on FDG', show=False)
        plt.savefig(os.path.join(out_dir, f"{group_ID}_qc_latenttime_fdg_{pcs}_{neighbours}_{count_level}_{paga_group}.png"),
                    dpi=300,
                    bbox_inches="tight")
        plt.close()

        #-----------------------------------------------------------------------------------
        # plot velocity arrows on both embeddings, coloured by cell cycle scores
        #-----------------------------------------------------------------------------------

        cc_score_tags = {"S.Score": "SScore", "G2M.Score": "G2MScore"}
        cc_embedding_stubs = {
            "umap_spliced":  f"{group_ID}_velocity_umap_{count_level}_{paga_group}",
            "draw_graph_fa": f"{group_ID}_velocity_fdg_{pcs}_{neighbours}_{count_level}_{paga_group}",
        }

        for score_col, score_tag in cc_score_tags.items():
            for basis, filename_stub in cc_embedding_stubs.items():
                scv.pl.velocity_embedding(
                    adata,
                    title=f"{score_col} - {basis} - dynamical - {count_level} - {group_ID}",
                    basis=basis,
                    color=score_col,
                    cmap='viridis',
                    arrow_length=3,
                    arrow_size=2,
                    dpi=120,
                    show=False
                )
                plt.savefig(os.path.join(out_dir, f"{filename_stub}_{score_tag}.png"),
                            dpi=300,
                            bbox_inches="tight")
                plt.close()

        #-----------------------------------------------------------------------------------
        # write out the velocity-processed object for downstream use
        #-----------------------------------------------------------------------------------

        adata.write_h5ad(
            os.path.join(out_dir, f"{group_ID}_velocity_dynamical_{pcs}_{neighbours}_{count_level}_{paga_group}.h5ad"),
            compression="gzip"
        )

        print(f"Finished: {group_ID}")

    except Exception as e:
        plt.close("all")  # make sure no partially-drawn figure is left open
        warnings.warn(f"{group_ID} — velocity run failed: {e}")
