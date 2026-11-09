#!/usr/bin/env python3
"""
Summarise PAGA connectivity across the 11.3 PC/neighbour sweep
================================================================

Loads every (pcs, neighbours) combination's .h5ad from
11.3-run_scanpy_paga_fdg_cell.py's output and reports the PAGA connectivity
value between Endocardium and Mesenchyme for each sample, so parameter
stability can be assessed at a glance rather than opening each saved PAGA
plot individually - some of which are skipped entirely when connectivity
falls below the plotting threshold (see 11.3), so a missing plot doesn't
mean the value can't be inspected here.

Also reports mean scVelo latent_time by cell type, from the single
velocity-processed object per sample (11.4-run_rnaVelocity_dynamical_fdg_
celltypes.py, fixed at pcs=30/neighbours=25) - latent_time is only computed
for that one combination, not across the full PAGA sweep above, so this is
a separate, sample-level summary rather than another row per combination.

Author: Melanie Smith
Email:  melanie.smith@adelaide.edu.au
Date:   2026-09-16

Usage:
    python summarise_paga_connectivity.py
"""

import os
import warnings

import pandas as pd
import scanpy as sc

#-----------------------------------------------------------------------------------
#   Params - match the sweep in 11.3-run_scanpy_paga_fdg_cell.py
#-----------------------------------------------------------------------------------

sample_names = ["e11_control", "e11_ko", "e12_control", "e12_ko"]
pc_range = [10, 15, 20, 30]
neighbour_range = [15, 25, 30, 50]
count_level = "spliced_unspliced"
paga_group = "celltypes"

base_dir = "/home/melanie-smith/workDir/sophieWiszniak/20251104_sophieWiszniak_ncc_pa/outDir/11-rna_velocity"

# Plotting threshold used in 11.3 - shown alongside each value for reference
plot_threshold = 0.01

#-----------------------------------------------------------------------------------
# Collect connectivity values across the sweep
#-----------------------------------------------------------------------------------

records = []

for sample_id in sample_names:
    for n_pcs in pc_range:
        for n_neighbours in neighbour_range:

            h5ad_path = os.path.join(
                base_dir, sample_id, "scanpy_output", paga_group,
                f"{sample_id}_fdg_pcs{n_pcs}_neighbours{n_neighbours}_counts_{count_level}_{paga_group}.h5ad"
            )

            connectivity = None

            if not os.path.exists(h5ad_path):
                warnings.warn(f"{sample_id} pcs={n_pcs} neighbours={n_neighbours}: file not found")
            else:
                try:
                    adata = sc.read_h5ad(h5ad_path)
                    conn_matrix = adata.uns["paga"]["connectivities"].toarray()
                    # off-diagonal entry - the single Endocardium/Mesenchyme edge
                    connectivity = conn_matrix[0, 1]
                except Exception as e:
                    warnings.warn(f"{sample_id} pcs={n_pcs} neighbours={n_neighbours}: {e}")

            records.append({
                "sample_id": sample_id,
                "n_pcs": n_pcs,
                "n_neighbours": n_neighbours,
                "connectivity": connectivity,
                "above_threshold": (connectivity is not None) and (connectivity >= plot_threshold),
            })

results = pd.DataFrame(records)

#-----------------------------------------------------------------------------------
# Print one pivoted table (pcs x neighbours) per sample
#-----------------------------------------------------------------------------------

pd.set_option("display.width", 120)

for sample_id in sample_names:
    print(f"\n=== {sample_id} — PAGA connectivity (Endocardium vs Mesenchyme) ===")
    sample_df = results[results["sample_id"] == sample_id]
    pivot = sample_df.pivot(index="n_pcs", columns="n_neighbours", values="connectivity")
    print(pivot.to_string(float_format=lambda x: f"{x:.4f}"))
    n_above = sample_df["above_threshold"].sum()
    print(f"{n_above} of {len(sample_df)} combinations >= threshold ({plot_threshold})")

#-----------------------------------------------------------------------------------
# Save the full long-format table for reference
#-----------------------------------------------------------------------------------

out_csv = os.path.join(base_dir, "paga_connectivity_summary.csv")
results.to_csv(out_csv, index=False)
print(f"\nSaved full summary to: {out_csv}")

#-----------------------------------------------------------------------------------
# Summarise latent_time by cell type, from the single velocity-processed object
# per sample (11.4-run_rnaVelocity_dynamical_fdg_celltypes.py; fixed at
# pcs=30, neighbours=25). Endocardium is expected to sit earlier (lower
# latent_time) than Mesenchyme if the inferred trajectory runs in the
# biologically expected EndMT direction.
#-----------------------------------------------------------------------------------

velocity_pcs = 30
velocity_neighbours = 25

latent_time_records = []

for sample_id in sample_names:

    velocity_h5ad_path = os.path.join(
        base_dir, sample_id, "velocity_output_dynamical", paga_group,
        f"{sample_id}_velocity_dynamical_pcs{velocity_pcs}_neighbours{velocity_neighbours}_{count_level}_{paga_group}.h5ad"
    )

    if not os.path.exists(velocity_h5ad_path):
        warnings.warn(f"{sample_id}: velocity h5ad not found at {velocity_h5ad_path}")
        continue

    try:
        adata = sc.read_h5ad(velocity_h5ad_path)
        mean_by_celltype = adata.obs.groupby("new_celltypes")["latent_time"].mean()

        latent_time_records.append({
            "sample_id": sample_id,
            "mean_latent_time_Endocardium": mean_by_celltype.get("Endocardium"),
            "mean_latent_time_Mesenchyme": mean_by_celltype.get("Mesenchyme"),
        })
    except Exception as e:
        warnings.warn(f"{sample_id}: failed to summarise latent_time: {e}")

latent_time_df = pd.DataFrame(latent_time_records)

if not latent_time_df.empty:
    latent_time_df["mesenchyme_minus_endocardium"] = (
        latent_time_df["mean_latent_time_Mesenchyme"] - latent_time_df["mean_latent_time_Endocardium"]
    )

print(f"\n=== Mean latent_time by cell type (pcs={velocity_pcs}, neighbours={velocity_neighbours}) ===")
print(latent_time_df.to_string(index=False, float_format=lambda x: f"{x:.4f}"))

latent_time_csv = os.path.join(base_dir, "latent_time_summary.csv")
latent_time_df.to_csv(latent_time_csv, index=False)
print(f"Saved latent_time summary to: {latent_time_csv}")
