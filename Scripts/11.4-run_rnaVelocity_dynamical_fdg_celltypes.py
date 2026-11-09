#!/usr/bin/env python3
"""
SCVELO
=================

Runs RNA velocity in dynamical mode on the PAGA/FDG-processed spliced/unspliced
Endocardium + Mesenchyme objects (from 11.3-run_scanpy_paga_fdg_cell.py), merging
in the raw Velocyto .loom spliced/unspliced/ambiguous counts.

Updated from run_rnaVelocity_dynamical_fdg_celltypes.py (archived) to:
  - point at the 11.3-generated PAGA/FDG .h5ad objects
    (outDir/11-rna_velocity/{sample_id}/scanpy_output/celltypes/...)
  - loop over all four samples, using a single fixed (pcs=30, neighbours=25)
    combination per sample, rather than a manually-edited sample_id/pcs/neighbours
  - drop the basis='umap' velocity plot - there is no gene-level UMAP in this
    object, only umap_spliced (confirmed against the actual .h5ad obsm keys)
  - fix the cluster column name used in the QC plot: full_dataset_clusters ->
    full_dataset_clusters_no_cc (confirmed actual column name)
  - add the missing scv.tl.latent_time(adata) call, required before the
    latent_time QC scatter plot (the archived template referenced 'latent_time'
    without ever computing it)
  - use n_jobs=4 for scv.tl.recover_dynamics() (running from the terminal, not
    interactively - see CLAUDE.md)
  - save the scv.pl.proportions() diagnostic plot to disk, matching how every
    other plot in this pipeline is handled, rather than only displaying it
  - export the finished velocity-processed adata object per sample as .h5ad -
    the archived template never saved this, only the plots
  - anchor scv.tl.latent_time()'s direction to a single, fixed root cell
    (the Endocardium cell most deeply embedded within the Endocardium
    population, by connectivity-graph density) rather than letting
    terminal_states() infer root/end purely from the velocity graph, which
    disagreed with the expected Endocardium -> Mesenchyme direction in most
    samples when left unanchored

Author: Melanie Smith
Email:  melanie.smith@adelaide.edu.au
Date:   2026-09-16

Usage:
    python 11.4-run_rnaVelocity_dynamical_fdg_celltypes.py

Notes:
    - Main libraries used: scvelo, scanpy, matplotlib
    - Python version: >= 3.9 recommended
"""

#-----------------------------------------------------------------------------------
#   Import core libraries
#-----------------------------------------------------------------------------------

import scvelo as scv
import scanpy as sc # single-cell analysis toolkit
import matplotlib.pyplot as plt
import numpy as np
import os
import warnings

#-----------------------------------------------------------------------------------
#   Params
#-----------------------------------------------------------------------------------

paga_group = "celltypes"
count_level = "spliced_unspliced"

sample_names = ["e11_control", "e11_ko", "e12_control", "e12_ko"]

# Fixed (pcs, neighbours) combination used for every sample, chosen as a single
# representative combination rather than repeating the expensive dynamical model
# fit across all 16 combinations swept in 11.3.
pcs_value = 30
neighbours_value = 25
pcs = f"pcs{pcs_value}"
neighbours = f"neighbours{neighbours_value}"

project_dir = "/home/melanie-smith/workDir/sophieWiszniak/20251104_sophieWiszniak_ncc_pa"

# Named vector mapping each sample to its Velocyto loom column prefix (same mapping
# used throughout the R/Python pipeline for this project)
loom_colnames = {
    "e11_control": "possorted_genome_bam_5Q1UF:",
    "e11_ko":      "possorted_genome_bam_LZRTJ:",
    "e12_control": "possorted_genome_bam_YIL7N:",
    "e12_ko":      "possorted_genome_bam_0YX73:",
}

#-----------------------------------------------------------------------------------
# Run RNA velocity (dynamical mode) for each sample
#-----------------------------------------------------------------------------------

for sample_id in sample_names:

    print(f"\n=== Processing: {sample_id} ===")

    try:
        # set the loom file address
        loom_file = os.path.join(
            project_dir, "outDir", "09-rna_velocity", sample_id, "velocyto_output",
            f"{loom_colnames[sample_id].rstrip(':')}.loom"
        )

        # set the .h5ad data input file address - the 11.3-generated PAGA/FDG object
        h5ad_file_input = os.path.join(
            project_dir, "outDir", "11-rna_velocity", sample_id, "scanpy_output", paga_group,
            f"{sample_id}_fdg_{pcs}_{neighbours}_counts_{count_level}_{paga_group}.h5ad"
        )

        # set the output directory
        out_dir = os.path.join(
            project_dir, "outDir", "11-rna_velocity", sample_id, "velocity_output_dynamical", paga_group
        )
        os.makedirs(out_dir, exist_ok=True)

        #-----------------------------------------------------------------------------------
        # Import the scv matrix and the intron/exon .loom file
        #-----------------------------------------------------------------------------------

        # import the normalised counts matrix from the scanpy_paga_fdg run (Seurat object with fdg added)
        adata = sc.read_h5ad(h5ad_file_input)

        # import the data from the 10x velocyto run (introns and exons)
        vlm = sc.read_loom(loom_file)

        #-----------------------------------------------------------------------------------
        # Inspect the .h5ad and .loom files -> do the cell barcodes match?
        #-----------------------------------------------------------------------------------

        print("h5ad example:", adata.obs_names[0])

        adata.obs_names = adata.obs_names.str.split("_").str[0]

        print(adata.obs_names[:5])

        print("loom example:", vlm.obs_names[0])

        # Remove sample prefix
        vlm.obs_names = vlm.obs_names.str.split(":").str[1]

        # Remove trailing x
        vlm.obs_names = vlm.obs_names.str.rstrip("x")

        # Add -1 suffix
        vlm.obs_names = vlm.obs_names + "-1"

        print(vlm.obs_names[:5])

        # check overlap
        common = adata.obs_names.intersection(vlm.obs_names)

        print("adata:", len(adata))
        print("loom:", len(vlm))
        print("common:", len(common))

        # subset the .loom file to match the adata object
        vlm = (vlm[common].copy())

        #-----------------------------------------------------------------------------------
        # merge the adata and ldata (vlm) objects
        #-----------------------------------------------------------------------------------

        # ensure cell order matches before merging
        vlm = vlm[adata.obs_names].copy()

        # merge the data
        adata = scv.utils.merge(adata, vlm)

        # safety check after merging
        print(adata.layers.keys())

        #-----------------------------------------------------------------------------------
        # let's take a look at the proportion of spliced and unspliced reads
        #-----------------------------------------------------------------------------------

        scv.pl.proportions(adata, show=False)
        plt.savefig(
            os.path.join(out_dir, f"{sample_id}_proportions_{pcs}_{neighbours}_{count_level}_{paga_group}.png"),
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
            title=f"umap_spliced projection - dynamical - {count_level} - {sample_id}",
            basis='umap_spliced',
            color='new_celltypes',
            arrow_length=3,
            arrow_size=2,
            dpi=120,
            show=False,
            legend_loc="bottom right"
        )
        plt.savefig(os.path.join(out_dir, f"{sample_id}_velocity_umap_{count_level}_{paga_group}.png"),
                    dpi=300,
                    bbox_inches="tight")
        plt.close()

        # QC plots
        scv.pl.scatter(adata, color='new_celltypes', basis='umap_spliced', title='Celltypes on UMAP', show=False)
        plt.savefig(os.path.join(out_dir, f"{sample_id}_qc_celltypes_umap_{count_level}_{paga_group}.png"),
                    dpi=300,
                    bbox_inches="tight")
        plt.close()

        scv.pl.scatter(adata, color='full_dataset_clusters_no_cc', basis='umap_spliced', title='Clusters on UMAP', show=False)
        plt.savefig(os.path.join(out_dir, f"{sample_id}_qc_clusters_umap_{count_level}_{paga_group}.png"),
                    dpi=300,
                    bbox_inches="tight")
        plt.close()

        scv.pl.scatter(adata, color='latent_time', basis='umap_spliced', cmap='gnuplot', title='Latent time on UMAP', show=False)
        plt.savefig(os.path.join(out_dir, f"{sample_id}_qc_latenttime_umap_{count_level}_{paga_group}.png"),
                    dpi=300,
                    bbox_inches="tight")
        plt.close()

        #-----------------------------------------------------------------------------------
        # plot velocity arrows on the FDG structure
        #-----------------------------------------------------------------------------------

        scv.pl.velocity_embedding(
            adata,
            title=f"PCs= {pcs}, Neighbours= {neighbours} - ({paga_group}) - {count_level} - dynamical - {sample_id}",
            basis='draw_graph_fa',
            color='new_celltypes',
            arrow_length=3,
            arrow_size=2,
            dpi=120,
            show=False,
            legend_loc="bottom right"
        )
        plt.savefig(os.path.join(out_dir, f"{sample_id}_velocity_fdg_{pcs}_{neighbours}_{count_level}_{paga_group}.png"),
                    dpi=300,
                    bbox_inches="tight")
        plt.close()

        # QC plots
        scv.pl.scatter(adata, color='new_celltypes', basis='draw_graph_fa', title='Celltypes on FDG', show=False)
        plt.savefig(os.path.join(out_dir, f"{sample_id}_qc_celltypes_fdg_{pcs}_{neighbours}_{count_level}_{paga_group}.png"),
                    dpi=300,
                    bbox_inches="tight")
        plt.close()

        scv.pl.scatter(adata, color='full_dataset_clusters_no_cc', basis='draw_graph_fa', title='Clusters on FDG', show=False)
        plt.savefig(os.path.join(out_dir, f"{sample_id}_qc_clusters_fdg_{pcs}_{neighbours}_{count_level}_{paga_group}.png"),
                    dpi=300,
                    bbox_inches="tight")
        plt.close()

        scv.pl.scatter(adata, color='latent_time', basis='draw_graph_fa', cmap='gnuplot', title='Latent time on FDG', show=False)
        plt.savefig(os.path.join(out_dir, f"{sample_id}_qc_latenttime_fdg_{pcs}_{neighbours}_{count_level}_{paga_group}.png"),
                    dpi=300,
                    bbox_inches="tight")
        plt.close()

        #-----------------------------------------------------------------------------------
        # write out the velocity-processed object for downstream use
        #-----------------------------------------------------------------------------------

        adata.write_h5ad(
            os.path.join(out_dir, f"{sample_id}_velocity_dynamical_{pcs}_{neighbours}_{count_level}_{paga_group}.h5ad"),
            compression="gzip"
        )

        print(f"Finished: {sample_id}")

    except Exception as e:
        plt.close("all")  # make sure no partially-drawn figure is left open
        warnings.warn(f"{sample_id} — velocity run failed: {e}")
