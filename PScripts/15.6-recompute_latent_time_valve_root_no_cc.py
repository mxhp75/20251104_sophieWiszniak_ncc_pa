#!/usr/bin/env python3
"""
Recompute latent_time with the root cell restricted to the Valve-scored
Endocardium population, per Sophie's request (project meeting follow-up,
2026-10-20).

15.x (no cc) counterpart of 14.7. Reads the already-computed 15.4 velocity
.h5ad and only redoes the root-cell choice + scv.tl.latent_time() - see
14.7 for the full rationale.

Valve labels: 14.5's endocardial_subtype scoring was only ever run on the
14.x (cc-retained) object - there is no 15.x counterpart. Reused directly
here rather than re-derived, since the 14.4/15.4 objects share an identical
14,631-cell population (confirmed by barcode set comparison) - Valve
identity is a property of the cell, not of whether cell-cycle genes were
later stripped out for this separate lineage's variance-partition workflow.

Overwrites:
  - all_samples_qc_latenttime_umap_spliced_unspliced_celltypes.png
  - all_samples_qc_latenttime_fdg_pcs30_neighbours25_spliced_unspliced_celltypes.png
  - all_samples_velocity_dynamical_pcs30_neighbours25_spliced_unspliced_celltypes.h5ad
    (same path as 15.4 - 15.5 reads this path downstream)

Author: Melanie Smith
Email:  melanie.smith@adelaide.edu.au
Date:   2026-10-20

Usage:
    python 15.6-recompute_latent_time_valve_root_no_cc.py
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

project_dir = "/home/melanie-smith/workDir/sophieWiszniak/20251104_sophieWiszniak_ncc_pa"
out_dir     = os.path.join(project_dir, "outDir", "15-rna_velocity")
velocity_out_dir = os.path.join(out_dir, "velocity_output_dynamical", paga_group)

h5ad_path = os.path.join(
    velocity_out_dir,
    f"all_samples_velocity_dynamical_{pcs}_{neighbours}_{count_level}_{paga_group}.h5ad"
)
# Valve labels only exist for the 14.x (cc-retained) object - reused here (see docstring)
valve_csv_path = os.path.join(
    project_dir, "outDir", "14-rna_velocity",
    "endocardial_mesenchymal_egfp", "all_samples_endocardial_subtype_assignments.csv"
)

print("\n=== Recomputing latent_time: Valve-restricted root (all_samples, no cell cycle) ===")
adata = sc.read_h5ad(h5ad_path)

#-----------------------------------------------------------------------------------
# Root cell: same density rule as 15.4, candidate pool narrowed to Valve-labelled cells
#-----------------------------------------------------------------------------------

subtype_df = pd.read_csv(valve_csv_path, index_col=0)
valve_barcodes = subtype_df.index[subtype_df["endocardial_subtype"] == "Valve"]
valve_in_object = adata.obs_names.isin(valve_barcodes)
print(f"Valve-labelled cells: {len(valve_barcodes)} total, {valve_in_object.sum()} present in this object")

root_mask = (valve_in_object & (adata.obs["new_celltypes"] == "Endocardium").values).astype(float)
root_density = adata.obsp["connectivities"].dot(root_mask)
root_density[root_mask == 0] = -np.inf
adata.uns["iroot"] = int(np.argmax(root_density))
print(f"New root cell: {adata.obs_names[adata.uns['iroot']]} (index {adata.uns['iroot']})")

scv.tl.latent_time(adata, root_key="iroot")

#-----------------------------------------------------------------------------------
# Overwrite the latent_time QC plots (UMAP + FDG)
#-----------------------------------------------------------------------------------

scv.pl.scatter(adata, color="latent_time", basis="umap_spliced", cmap="gnuplot",
                title="Latent time on UMAP - Valve root (no cc) - all_samples", show=False)
plt.savefig(os.path.join(velocity_out_dir, f"all_samples_qc_latenttime_umap_{count_level}_{paga_group}.png"),
            dpi=300, bbox_inches="tight")
plt.close()

scv.pl.scatter(adata, color="latent_time", basis="draw_graph_fa", cmap="gnuplot",
                title="Latent time on FDG - Valve root (no cc) - all_samples", show=False)
plt.savefig(os.path.join(velocity_out_dir, f"all_samples_qc_latenttime_fdg_{pcs}_{neighbours}_{count_level}_{paga_group}.png"),
            dpi=300, bbox_inches="tight")
plt.close()

#-----------------------------------------------------------------------------------
# Overwrite the velocity object (same path 15.5 reads) so downstream plots pick up
# the new latent_time
#-----------------------------------------------------------------------------------

adata.write_h5ad(h5ad_path, compression="gzip")

print("Finished: all_samples (no cell cycle)")
