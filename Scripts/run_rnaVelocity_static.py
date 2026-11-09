#!/usr/bin/env python3
"""
SCVELO
=================

This script runs RNA-Velocity
- Prior to running RNA-Velocity the .loom file containing the intron/exon needs to be created

Author: Melanie Smith
Email:  melanie.smith@adelaide.edu.au
Date:   2026-02-17

Usage:
    python this_script.py [options] [arguments]

Example:
    python this_script.py --input data.csv --output results/

Notes:
    - Main libraries used: ...
    - Python version: >= 3.9 recommended
"""
from operator import index

#-----------------------------------------------------------------------------------
#   Global configs
#-----------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------
#   Import core libraries
#-----------------------------------------------------------------------------------
import scvelo as scv
import velocyto as vcy
import pandas as pd
import scanpy as sc # single-cell analysis toolkit
import matplotlib.pyplot as plt
import os

#-----------------------------------------------------------------------------------
#   Params
#-----------------------------------------------------------------------------------

# set the sample ID
sample_id = "e11_control"
#sample_id = "e11_ko"
#sample_id = "e12_control"
#sample_id = "e12_ko"

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

# set the .h5ad file address
five_hd_file = os.path.join(
    project_dir,
    "outDir",
    "09-rna_velocity",
    sample_id,
    f"{sample_id}_fdg.h5ad")

# set the output directory
out_dir = os.path.join(
    project_dir,
    "outDir",
    "09-rna_velocity",
    sample_id,
    "velocity_output"
)
os.makedirs(out_dir, exist_ok=True)

#-----------------------------------------------------------------------------------
# Import the scv matrix and the intron/exon .loom file
#-----------------------------------------------------------------------------------

# import the normalised counts matrix from the scanpy_paga_fdg run (old Seurat object with dfg added)
adata = sc.read_h5ad(five_hd_file)

# import the data from the 10x velocyto run (introns and exons)
vlm = sc.read_loom(loom_file)

#-----------------------------------------------------------------------------------
# merge the adata and ldata (vlm) objects
#-----------------------------------------------------------------------------------

adata = scv.utils.merge(adata, vlm)

#-----------------------------------------------------------------------------------
# let's take a look at the proportion of spliced and unspliced reads
#-----------------------------------------------------------------------------------
scv.pl.proportions(adata)

#-----------------------------------------------------------------------------------
# now we need to normalise the spliced data (Seurat object is already normalised - should detect this but need to check)
#-----------------------------------------------------------------------------------
scv.pp.filter_and_normalize(adata)

#-----------------------------------------------------------------------------------
# calculate the moments (i.e. means and uncentered variances) computed among nearest neighbors in PCA space
#-----------------------------------------------------------------------------------
scv.pp.moments(adata)

#-----------------------------------------------------------------------------------
# compute the velocities
#-----------------------------------------------------------------------------------
scv.tl.velocity(adata)

#-----------------------------------------------------------------------------------
# compute the velocity graph
#-----------------------------------------------------------------------------------
scv.tl.velocity_graph(adata)

#-----------------------------------------------------------------------------------
# plot each individual cell as one velocity arrow using the umap projections
#-----------------------------------------------------------------------------------
scv.pl.velocity_embedding(
    adata,
    title=f"umap projection - steady state - {sample_id}",
    basis='umap',
    color='new_celltypes',
    arrow_length=3,
    arrow_size=2,
    dpi=120,
    show=False,
    legend_loc="bottom right"
)
plt.savefig(os.path.join(out_dir, f"{sample_id}_velocity_umap.png"),
            dpi=300,
            bbox_inches="tight")
plt.close()

# plot each individual cell as one velocity arrow using the fdg projections
scv.pl.velocity_embedding(
    adata,
    title=f"FDG projection - steady state - {sample_id}",
    basis='draw_graph_fa',
    color='new_celltypes',
    arrow_length=3,
    arrow_size=2,
    dpi=120,
    show=False,
    legend_loc="bottom right"
)
plt.savefig(os.path.join(out_dir, f"{sample_id}_velocity_fdg.png"),
            dpi=300,
            bbox_inches="tight")
plt.close()

# plot the embedding on a grid -> better visualisation
scv.pl.velocity_embedding_grid(adata,
                               basis='umap',
                               color='new_celltypes')

# plot the umap projections embedding as a stream
scv.pl.velocity_embedding_stream(adata,
                                 title=f"umap projection - steady state - {sample_id}",
                                 basis='umap',
                                 color='new_celltypes',
                                 show=False,
                                 legend_loc="bottom right"
                                 )
plt.savefig(os.path.join(out_dir, f"{sample_id}_velocity_stream_umap.png"),
            dpi=300,
            bbox_inches="tight")
plt.close()

# plot the FDG projections embedding as a stream
scv.pl.velocity_embedding_stream(adata,
                                 title=f"FDG projection - steady state - {sample_id}",
                                 basis='draw_graph_fa',
                                 color='new_celltypes',
                                 show=False,
                                 legend_loc="bottom right"
                                 )
plt.savefig(os.path.join(out_dir, f"{sample_id}_velocity_stream_fdg.png"),
            dpi=300,
            bbox_inches="tight")
plt.close()


 #-----------------------------------------------------------------------------------
 # examine the phase portraits of marker genes
 #-----------------------------------------------------------------------------------
scv.pl.velocity(adata,
                ['Postn',  'Cthrc1'],
                basis='draw_graph_fa',
                color='new_celltypes',
                ncols=1,
                nrows=2,
                show=False,
                dpi=300)
plt.savefig(os.path.join(out_dir, f"{sample_id}_phase_portrait_fdg_mesenchyme.png"),
            dpi=300,
            bbox_inches="tight")
plt.close()

scv.pl.velocity(adata,
                ['Postn',  'Cthrc1'],
                basis='umap',
                color='new_celltypes',
                ncols=1,
                nrows=2,
                show=False,
                dpi=300)
plt.savefig(os.path.join(out_dir, f"{sample_id}_phase_portrait_umap_mesenchyme.png"),
            dpi=300,
            bbox_inches="tight")
plt.close()