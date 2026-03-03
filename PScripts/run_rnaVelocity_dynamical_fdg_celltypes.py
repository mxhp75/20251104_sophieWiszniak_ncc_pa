#!/usr/bin/env python3
"""
SCVELO
=================

This script runs RNA-Velocity in dynamical mode

Author: Melanie Smith
Email:  melanie.smith@adelaide.edu.au
Date:   2026-02-27

Usage:
    python this_script.py [options] [arguments]

Example:
    python this_script.py --input data.csv --output results/

Notes:
    - Main libraries used: ...
    - Python version: >= 3.9 recommended
"""

#-----------------------------------------------------------------------------------
#   Import core libraries
#-----------------------------------------------------------------------------------
import scvelo as scv
import scanpy as sc # single-cell analysis toolkit
import matplotlib.pyplot as plt
import os

#-----------------------------------------------------------------------------------
#   Params
#-----------------------------------------------------------------------------------

# set PAGA group
paga_group = "celltypes"
# paga_group = "clusters"

# set the sample ID
sample_id = "e11_control"
#sample_id = "e11_ko"
#sample_id = "e12_control"
#sample_id = "e12_ko"

# set the number of PCs
#pcs="pcs20"
pcs="pcs30"

# set the number of neighbours
#neighbours="neighbours15"
neighbours="neighbours25"
#neighbours="neighbours50"

# set the main project directory
project_dir = f"/home/melanie-smith/workDir/sophieWiszniak/20251104_sophieWiszniak_ncc_pa"

# set the loom file address
loom_file = os.path.join(
    project_dir,
    "outDir",
    "09-rna_velocity",
    sample_id,
    "velocyto_output",
    "possorted_genome_bam_5Q1UF.loom"
)

# set the .h5ad data input file address
h5ad_file_input = os.path.join(
    project_dir,
    "outDir",
    "09-rna_velocity",
    sample_id,
    "scanpy_output",
    paga_group,
    f"{sample_id}_fdg_{pcs}_{neighbours}_{paga_group}.h5ad")

# set the output directory
out_dir = os.path.join(
    project_dir,
    "outDir",
    "09-rna_velocity",
    sample_id,
    "velocity_output_dynamical",
    paga_group
)
os.makedirs(out_dir, exist_ok=True)

#-----------------------------------------------------------------------------------
# Import the scv matrix and the intron/exon .loom file
#-----------------------------------------------------------------------------------

# import the normalised counts matrix from the scanpy_paga_fdg run (old Seurat object with fdg added)
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
scv.pl.proportions(adata)

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

scv.tl.recover_dynamics(adata, n_jobs=1, show_progress_bar=False)          # this can take time
scv.tl.velocity(adata, mode='dynamical')

#-----------------------------------------------------------------------------------
# compute the velocity graph
#-----------------------------------------------------------------------------------

scv.tl.velocity_graph(adata)

#-----------------------------------------------------------------------------------
# plot each individual cell as one velocity arrow using the umap projections
#-----------------------------------------------------------------------------------
scv.pl.velocity_embedding(
    adata,
    title=f"umap projection - dynamical - {sample_id}",
    basis='umap',
    color='new_celltypes',
    arrow_length=3,
    arrow_size=2,
    dpi=120,
    show=False,
    legend_loc="bottom right"
)
plt.savefig(os.path.join(out_dir, f"{sample_id}_velocity_umap_{pcs}_{neighbours}_{paga_group}.png"),
            dpi=300,
            bbox_inches="tight")
plt.close()

scv.pl.velocity_embedding(
    adata,
    title=f"PCs= {pcs}, Neighbours= {neighbours} - ({paga_group}) - dynamical - {sample_id}",
    basis='draw_graph_fa',
    color='new_celltypes',
    arrow_length=3,
    arrow_size=2,
    dpi=120,
    show=False,
    legend_loc="bottom right"
)
plt.savefig(os.path.join(out_dir, f"{sample_id}_velocity_fdg_{pcs}_{neighbours}_{paga_group}.png"),
            dpi=300,
            bbox_inches="tight")
plt.close()
