#!/usr/bin/env python3
"""
SCVELO
=================

This script imports the 4 x andata objects, one per sample, and explores the velocity and latent tine plots

Author: Melanie Smith
Email:  melanie.smith@adelaide.edu.au
Date:   2026-03-19

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
from pathlib import Path



#-----------------------------------------------------------------------------------
#   adata files to import
#-----------------------------------------------------------------------------------

# set the main project directory
project_dir = f"/home/melanie-smith/workDir/sophieWiszniak/20251104_sophieWiszniak_ncc_pa"

adata_e11_control_file = os.path.join(
    project_dir,
    "outDir", "09-rna_velocity", "e11_control", "velocity_output_dynamical", "clusters",
    f"e11_control_velocity_fdg_pcs20_neighbours50_spliced_unspliced_clusters.h5ad"
    )

adata_e11_ko_file = os.path.join(
    project_dir,
    "outDir", "09-rna_velocity", "e11_ko", "velocity_output_dynamical", "clusters",
    f"e11_ko_velocity_fdg_pcs20_neighbours50_spliced_unspliced_clusters.h5ad"
    )

adata_e12_control_file = os.path.join(
    project_dir,
    "outDir", "09-rna_velocity", "e12_control", "velocity_output_dynamical", "clusters",
    f"e12_control_velocity_fdg_pcs30_neighbours15_spliced_unspliced_clusters.h5ad"
    )

adata_e12_ko_file = os.path.join(
    project_dir,
    "outDir", "09-rna_velocity", "e12_ko", "velocity_output_dynamical", "clusters",
    f"e12_ko_velocity_fdg_pcs20_neighbours50_spliced_unspliced_clusters.h5ad"
    )

#-----------------------------------------------------------------------------------
# Import the sample-specific .h5ad files
#-----------------------------------------------------------------------------------

# import the normalised counts matrix from the scanpy_paga_fdg run (old Seurat object with fdg added)
adata_e11_control = sc.read_h5ad(adata_e11_control_file)
adata_e11_ko = sc.read_h5ad(adata_e11_ko_file)
adata_e12_control = sc.read_h5ad(adata_e12_control_file)
adata_e12_ko = sc.read_h5ad(adata_e12_ko_file)

#-----------------------------------------------------------------------------------
# Generate the QC plots and save
#-----------------------------------------------------------------------------------

# Heatmap top 10
## make a list of adata objects
base_dir = Path(project_dir)

adatas = {
    "e11_control": adata_e11_control,
    "e11_ko": adata_e11_ko,
    "e12_control": adata_e12_control,
    "e12_ko": adata_e12_ko
}

# heatmap of top genes against latent time
for sample_id, adata in adatas.items():

    out_dir = base_dir / "outDir" / "09-rna_velocity" / sample_id / "velocity_output_dynamical" / "clusters"
    out_dir.mkdir(parents=True, exist_ok=True)

    adata = adata.copy()

    top_genes = (
        adata.var['fit_likelihood']
        .dropna()
        .sort_values(ascending=False)
        .index[:300]
    )

    cg = scv.pl.heatmap(
        adata,
        var_names=top_genes,
        sortby='latent_time',
        col_color='new_celltypes',
        n_convolve=100,
        show=False
    )

    out_file = out_dir / f"{sample_id}_heatmap_300.png"
    cg.savefig(out_file, dpi=300)

    plt.close('all')

    print(f"Saved: {out_file}")

# Velocity plot using the cluster based paga seed but coloured by the celltypes

for sample_id, adata in adatas.items():

    project_dir = Path("/home/melanie-smith/workDir/sophieWiszniak/20251104_sophieWiszniak_ncc_pa")
    out_dir = project_dir / "outDir" / "09-rna_velocity" / sample_id / "velocity_output_dynamical" / "clusters"
    out_dir.mkdir(parents=True, exist_ok=True)

    adata = adata.copy()

    vp = scv.pl.velocity_embedding(
        adata,
        title=f"Velocity - dynamical - {sample_id}",
        basis='draw_graph_fa',
        color='new_celltypes',
        arrow_length=3,
        arrow_size=2,
        dpi=120,
        show=False,
        legend_loc="bottom right"
    )

    out_file = out_dir / f"{sample_id}_velocity_cellType_colour.png"
    vp.get_figure().savefig(out_file, dpi=300, bbox_inches='tight')
    plt.close('all')
    print(f"Saved: {out_file}")

# Phase diagrams in a loop
for sample_id, adata in adatas.items():

    out_dir = project_dir / "outDir" / "09-rna_velocity" / sample_id / "velocity_output_dynamical" / "clusters"
    out_dir.mkdir(parents=True, exist_ok=True)

    adata = adata.copy()

    top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index

    def save_scatter(ax, out_file):
        fig = ax[0].get_figure() if isinstance(ax, list) else ax.get_figure()
        fig.savefig(out_file, dpi=300, bbox_inches='tight')
        plt.close('all')
        print(f"Saved: {out_file}")

    # phase diagram of top likelihood genes coloured by clusters
    ax = scv.pl.scatter(adata,
                        basis=top_genes[:15],
                        ncols=5,
                        frameon=True,
                        color='full_dataset_clusters',
                        dpi=90,
                        show=False
                        )
    save_scatter(ax, out_dir / f"{sample_id}_phase_diagram_top15_clusters.png")

    # phase diagram of top likelihood genes coloured by cell type
    ax = scv.pl.scatter(adata,
                        basis=top_genes[:15],
                        ncols=5,
                        frameon=True,
                        color='new_celltypes',
                        dpi=90,
                        show=False
                        )
    save_scatter(ax, out_dir / f"{sample_id}_phase_diagram_top15_celltype.png")

    # mesenchymal genes
    var_names_mesenchymal = ['Postn', 'Vcan', 'Dbi', 'Cthrc1', 'Papss2']
    ax = scv.pl.scatter(adata,
                        var_names_mesenchymal,
                        color='full_dataset_clusters',
                        dpi=120,
                        show=False,
                        legend_loc='none'
                        )
    save_scatter(ax, out_dir / f"{sample_id}_phase_diagram_mesenchymal_clusters.png")

    ax = scv.pl.scatter(adata,
                        var_names_mesenchymal,
                        color='new_celltypes',
                        legend_loc='none',
                        show=False
                        )
    save_scatter(ax, out_dir / f"{sample_id}_phase_diagram_mesenchymal_celltype.png")

    ax = scv.pl.scatter(adata,
                        x='latent_time',
                        y=var_names_mesenchymal,
                        color='new_celltypes',
                        dpi=120,
                        frameon=False,
                        show=False,
                        legend_loc='none'
                        )
    save_scatter(ax, out_dir / f"{sample_id}_scatter_latentTime_expression_mesenchymal.png")

    # epicardial genes
    var_names_epicardial = ['Wt1', 'Upk3b', 'Apela', 'Upk1b', 'Aldh1a2']
    ax = scv.pl.scatter(adata,
                        var_names_epicardial,
                        color='full_dataset_clusters',
                        dpi=120,
                        show=False,
                        legend_loc='none'
                        )
    save_scatter(ax, out_dir / f"{sample_id}_phase_diagram_epicardial_clusters.png")

    ax = scv.pl.scatter(adata,
                        var_names_epicardial,
                        color='new_celltypes',
                        dpi=120,
                        legend_loc='none',
                        show=False
                        )
    save_scatter(ax, out_dir / f"{sample_id}_phase_diagram_epicardial_celltype.png")

    ax = scv.pl.scatter(adata,
                        x='latent_time',
                        y=var_names_epicardial,
                        color='new_celltypes',
                        dpi=120,
                        frameon=False,
                        show=False,
                        legend_loc='none'
                        )
    save_scatter(ax, out_dir / f"{sample_id}_scatter_latentTime_expression_epicardial.png")

    # endocardial genes
    var_names_endocardial = ['Ecscr', 'Cdh5', 'Tie1', 'Rasip1', 'Plxnd1']
    ax = scv.pl.scatter(adata,
                        var_names_endocardial,
                        color='full_dataset_clusters',
                        dpi=120,
                        frameon=False,
                        show=False,
                        legend_loc='none'
                        )
    save_scatter(ax, out_dir / f"{sample_id}_phase_diagram_endocardial_clusters.png")

    ax = scv.pl.scatter(adata,
                        var_names_endocardial,
                        color='new_celltypes',
                        dpi=120,
                        legend_loc='none',
                        frameon=False,
                        show=False
                        )
    save_scatter(ax, out_dir / f"{sample_id}_phase_diagram_endocardial_celltype.png")

    ax = scv.pl.scatter(adata,
                        x='latent_time',
                        y=var_names_endocardial,
                        color='new_celltypes',
                        dpi=120,
                        frameon=False,
                        show=False,
                        legend_loc='none'
                        )
    save_scatter(ax, out_dir / f"{sample_id}_scatter_latentTime_expression_endocardial.png")

    # VSMC genes
    var_names_vsmc = ['Cxcl12', 'Iigp1', 'Eln', 'Fbln5', 'Rgs5']
    ax = scv.pl.scatter(adata,
                        var_names_vsmc,
                        color='full_dataset_clusters',
                        dpi=120,
                        frameon=False,
                        show=False,
                        legend_loc='none'
                        )
    save_scatter(ax, out_dir / f"{sample_id}_phase_diagram_vsmc_clusters.png")

    ax = scv.pl.scatter(adata,
                        var_names_vsmc,
                        color='new_celltypes',
                        dpi=120,
                        legend_loc='none',
                        frameon=False,
                        show=False
                        )
    save_scatter(ax, out_dir / f"{sample_id}_phase_diagram_vsmc_celltype.png")

    ax = scv.pl.scatter(adata,
                        x='latent_time',
                        y=var_names_vsmc,
                        color='new_celltypes',
                        dpi=120,
                        frameon=False,
                        show=False,
                        legend_loc='none'
                        )
    save_scatter(ax, out_dir / f"{sample_id}_scatter_latentTime_expression_vsmc.png")

    # myocardial genes
    var_names_myocardial = ['Chchd10', 'Sh3bgr', 'Ttn', 'Pln', 'Myh7']
    ax = scv.pl.scatter(adata,
                        var_names_myocardial,
                        color='full_dataset_clusters',
                        dpi=120,
                        frameon=False,
                        show=False,
                        legend_loc='none'
                        )
    save_scatter(ax, out_dir / f"{sample_id}_phase_diagram_myocardial_clusters.png")

    ax = scv.pl.scatter(adata,
                        var_names_myocardial,
                        color='new_celltypes',
                        dpi=120,
                        legend_loc='none',
                        frameon=False,
                        show=False
                        )
    save_scatter(ax, out_dir / f"{sample_id}_phase_diagram_myocardial_celltype.png")

    ax = scv.pl.scatter(adata,
                        x='latent_time',
                        y=var_names_myocardial,
                        color='new_celltypes',
                        dpi=120,
                        frameon=False,
                        show=False,
                        legend_loc='none'
                        )
    save_scatter(ax, out_dir / f"{sample_id}_scatter_latentTime_expression_myocardial.png")


# Phase diagrams of genes Sophie is interested in
scv.pl.scatter(
    adata_e12_ko,
    basis=['Plxna2', 'Itga4', 'Irx5', 'Tgfbr3', 'Apln', 'Sema3d', 'Sema3c', 'Vegfa'],
    ncols=5,
    frameon=True,
    color='new_celltypes',
    dpi=90,
    show=False
)

plt.suptitle("E12 ko", fontsize=16)
plt.show()
