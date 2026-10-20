#!/usr/bin/env python3
"""
Recompute latent_time with the root cell restricted to the Valve-scored
Endocardium population, per Sophie's request (mirrors 14.7/15.6, applied
here to the separately-fitted genotype-split models used for slide 19 of
the lab meeting deck).

Not a new dynamical-model fit - reads the already-computed 13.4 velocity
.h5ad per group (recover_dynamics/velocity/velocity_graph already done
there) and only redoes the root-cell choice + scv.tl.latent_time(), which
is the sole step downstream of the root. Mirrors 13.4's root-selection
algorithm exactly (the Endocardium cell with the highest connectivity-graph
density of other Endocardium cells), just with the candidate pool narrowed
to Valve-labelled cells instead of all Endocardium cells.

Valve labels come from 13.6's own per-group endocardial_subtype assignment
(each genotype group has its own independently-fitted subtype labels, not
shared with the 14.x/15.x merged-model labels used in 14.7/15.6).

Overwrites, per group:
  - {group}_qc_latenttime_umap_spliced_unspliced_celltypes.png
  - {group}_qc_latenttime_fdg_pcs30_neighbours25_spliced_unspliced_celltypes.png
  - {group}_velocity_dynamical_pcs30_neighbours25_spliced_unspliced_celltypes.h5ad

Author: Melanie Smith
Email:  melanie.smith@adelaide.edu.au
Date:   2026-10-20

Usage:
    python 13.8-recompute_latent_time_valve_root.py
"""

#-----------------------------------------------------------------------------------
#   Import core libraries
#-----------------------------------------------------------------------------------

import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
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

group_names = ["control", "ko"]

project_dir = "/home/melanie-smith/workDir/sophieWiszniak/20251104_sophieWiszniak_ncc_pa"

for group_ID in group_names:

    print(f"\n=== Recomputing latent_time: Valve-restricted root ({group_ID}) ===")

    out_dir = os.path.join(
        project_dir, "outDir", "13-rna_velocity", group_ID, "velocity_output_dynamical", paga_group
    )
    h5ad_path = os.path.join(
        out_dir, f"{group_ID}_velocity_dynamical_{pcs}_{neighbours}_{count_level}_{paga_group}.h5ad"
    )
    valve_csv_path = os.path.join(
        project_dir, "outDir", "13-rna_velocity", group_ID, "endocardial_subtypes",
        f"{group_ID}_endocardial_subtype_assignments.csv"
    )

    adata = sc.read_h5ad(h5ad_path)

    #-----------------------------------------------------------------------------------
    # Root cell: same density rule as 13.4, candidate pool narrowed to this group's
    # own Valve-labelled cells
    #-----------------------------------------------------------------------------------

    subtype_df = pd.read_csv(valve_csv_path, index_col=0)
    valve_barcodes = subtype_df.index[subtype_df["endocardial_subtype"] == "Valve"]
    valve_in_object = adata.obs_names.isin(valve_barcodes)
    print(f"Valve-labelled cells ({group_ID}): {len(valve_barcodes)} total, {valve_in_object.sum()} present in this object")

    root_mask = (valve_in_object & (adata.obs["new_celltypes"] == "Endocardium").values).astype(float)
    root_density = adata.obsp["connectivities"].dot(root_mask)
    root_density[root_mask == 0] = -np.inf
    adata.uns["iroot"] = int(np.argmax(root_density))
    print(f"New root cell ({group_ID}): {adata.obs_names[adata.uns['iroot']]} (index {adata.uns['iroot']})")

    scv.tl.latent_time(adata, root_key="iroot")

    #-----------------------------------------------------------------------------------
    # Overwrite the latent_time QC plots (UMAP + FDG)
    #-----------------------------------------------------------------------------------

    scv.pl.scatter(adata, color="latent_time", basis="umap_spliced", cmap="gnuplot",
                    title=f"Latent time on UMAP - Valve root - {group_ID}", show=False)
    plt.savefig(os.path.join(out_dir, f"{group_ID}_qc_latenttime_umap_{count_level}_{paga_group}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()

    scv.pl.scatter(adata, color="latent_time", basis="draw_graph_fa", cmap="gnuplot",
                    title=f"Latent time on FDG - Valve root - {group_ID}", show=False)
    plt.savefig(os.path.join(out_dir, f"{group_ID}_qc_latenttime_fdg_{pcs}_{neighbours}_{count_level}_{paga_group}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()

    #-----------------------------------------------------------------------------------
    # Overwrite the velocity object (same path 13.6/other downstream scripts read)
    #-----------------------------------------------------------------------------------

    adata.write_h5ad(h5ad_path, compression="gzip")

    print(f"Finished: {group_ID}")
