#!/usr/bin/env python3
"""

=================

This script takes the spliced/unspliced Endocardium + Mesenchyme object
(re-embedded, converted to .h5ad by 14.2-loom_to_seurat_endo_mesen.Rmd) and
generates the PAGA graph (coarse-grained connectivity map between cell type
clusters), then uses this information to generate the FDG projections.

Mirrors 12.3/13.3-run_scanpy_paga_fdg_cell.py, updated to:
  - point at the 14.2-generated spliced/unspliced .h5ad object
    (outDir/14-rna_velocity/all_samples_spliced_endo_mesen_subset.h5ad) -
    the all-four-samples-merged, EGFP-lineage-tracing counterpart of
    12.2/13.2's genotype-pooled objects
  - no outer loop over groups/samples - this is a single merged object, so
    the (n_pcs, n_neighbours) sweep runs directly on it

Author: Melanie Smith
Email:  melanie.smith@adelaide.edu.au
Date:   2026-10-07

Usage:
    python 14.3-run_scanpy_paga_fdg_cell.py

Notes:
    - Main libraries used: scanpy, matplotlib
    - Python version: >= 3.9 recommended
"""

#-----------------------------------------------------------------------------------
#   Import core libraries
#-----------------------------------------------------------------------------------

import scanpy as sc # single-cell analysis toolkit
import matplotlib.pyplot as plt
import os
import warnings

#-----------------------------------------------------------------------------------
#   Params
#-----------------------------------------------------------------------------------

fPCA = 50  # no. of PCs to compute once - covers all tested n_pcs values below

project_dir = "/home/melanie-smith/workDir/sophieWiszniak/20251104_sophieWiszniak_ncc_pa"
out_dir     = os.path.join(project_dir, "outDir", "14-rna_velocity")

h5ad_input_path = os.path.join(out_dir, "all_samples_spliced_endo_mesen_subset.h5ad")

# Label used only in output filenames, matching the existing project convention
# (see CLAUDE.md: "..._counts_{count_level}_celltypes..."). Fixed here since the
# 14.2 object is the spliced/unspliced subset - not a variable choice.
count_level = "spliced_unspliced"
paga_group  = "celltypes"

pc_range = [10, 15, 20, 30]         # PCs considered - sliced from the single fPCA-component PCA
neighbour_range = [15, 25, 30, 50]

#-----------------------------------------------------------------------------------
# Create embeddings for dynamic velocity modelling - PAGA informed Force Directed Graph
# Single merged object (all four samples) - no outer loop over groups/samples.
# Save output for all permutations.
#-----------------------------------------------------------------------------------

print(f"\n=== Processing: all_samples ===")

adata = sc.read_h5ad(h5ad_input_path)
print(adata)
print(adata.layers.keys())       # at least 'counts_RNA'
print(adata.obsm.keys())         # X_pca_spliced, X_umap_spliced

# new_celltypes still carries the full cell-type category list from before the
# 14.1 Endocardium/Mesenchyme subsetting (e.g. Myocardium, VSMC, Epicardium),
# even though only a couple of those categories have cells in this object. PAGA
# sizes its group-position array off the declared categories, so leftover unused
# levels cause an out-of-bounds index in sc.pl.paga() - drop them first.
adata.obs["new_celltypes"] = adata.obs["new_celltypes"].cat.remove_unused_categories()

# Compute PCA once with enough components (50) to cover all tested n_pcs values
sc.tl.pca(adata, svd_solver='arpack', n_comps=fPCA)

# Create output directory if it doesn't exist
save_dir = os.path.join(out_dir, "scanpy_output", paga_group)
os.makedirs(save_dir, exist_ok=True)

n_failed = 0
n_paga_plot_skipped = 0

for n_pcs in pc_range:
    for n_neighbours in neighbour_range:

        # Each (n_pcs, n_neighbours) combination is wrapped so that a single
        # failure doesn't halt the sweep for the rest of the combinations.
        try:
            adata_test = adata.copy()

            sc.pp.neighbors(adata_test, n_pcs=n_pcs, n_neighbors=n_neighbours)

            # Compute PAGA connectivities - kept regardless of whether the
            # diagnostic plot below succeeds, since the connectivity value
            # itself is a useful result in its own right.
            sc.tl.paga(adata_test, groups="new_celltypes")

            # Plot & SAVE PAGA - isolated in its own try/except. With only two
            # broad cell-type groups, the connectivity between them can be at
            # or below the plotting threshold (a real "no significant
            # connectivity" result rather than an error), which crashes
            # sc.pl.paga() in this scanpy version. That failure should only
            # cost us the PAGA diagnostic plot, not the FDG/export below -
            # FDG is now seeded from PCA coordinates rather than PAGA
            # positions, so it no longer depends on this plot succeeding.
            try:
                sc.pl.paga(
                    adata_test,
                    threshold=0.01,
                    color="new_celltypes",
                    node_size_scale=4,
                    fontsize=5,
                    title=f"PAGA – all_samples – pcs={n_pcs}, neighbours={n_neighbours}, counts={count_level}",
                    show=False
                )

                plt.savefig(
                    os.path.join(save_dir, f"paga_pcs{n_pcs}_neighbours{n_neighbours}_counts_{count_level}_celltypes.png"),
                    bbox_inches="tight",
                    dpi=150
                )
                plt.close()  # Clean up this figure immediately

            except Exception as e:
                n_paga_plot_skipped += 1
                plt.close("all")
                warnings.warn(
                    f"all_samples — pcs={n_pcs}, neighbours={n_neighbours}: "
                    f"PAGA plot skipped (likely no significant connectivity): {e}"
                )

            # Compute and plot FDG - seeded from the existing spliced PCA
            # coordinates rather than PAGA positions. With only two broad
            # cell-type groups, PAGA-seeding offers little benefit over a
            # PCA-based start.
            sc.tl.draw_graph(adata_test, init_pos="X_pca")
            sc.pl.embedding(
                adata_test,
                basis='draw_graph_fa',  # explicit basis (no X_ prefix needed)
                color="new_celltypes",
                title=f"FDG – all_samples – pcs={n_pcs}, neighbours={n_neighbours}, counts={count_level}",
                legend_loc="right margin",
                legend_fontsize=8,
                size=18,
                frameon=False,
                show=False
            )

            plt.savefig(
                os.path.join(save_dir, f"fdg_pcs{n_pcs}_neighbours{n_neighbours}_counts_{count_level}_celltypes.png"),
                bbox_inches="tight",
                dpi=150
            )
            plt.close()

            # write out the adata object for use in the RNA Velocity workflow
            adata_test.write_h5ad(
                os.path.join(save_dir, f"all_samples_fdg_pcs{n_pcs}_neighbours{n_neighbours}_counts_{count_level}_celltypes.h5ad"),
                compression="gzip"
            )

        except Exception as e:
            n_failed += 1
            plt.close("all")  # make sure no partially-drawn figure is left open
            warnings.warn(
                f"all_samples — pcs={n_pcs}, neighbours={n_neighbours} failed: {e}"
            )

n_total = len(pc_range) * len(neighbour_range)
print(
    f"Finished: all_samples ({n_failed} of {n_total} combinations failed, "
    f"{n_paga_plot_skipped} of {n_total} PAGA plots skipped)"
)
