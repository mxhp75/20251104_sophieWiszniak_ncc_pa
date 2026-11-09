#!/usr/bin/env python3

"""
=================

This script runs RNA-Velocity in dynamical mode

Author: Melanie Smith
Email:  melanie.smith@adelaide.edu.au
Date:   2026-02-27

=================
"""

#-----------------------------------------------------------------------------------
#   Import core libraries
#-----------------------------------------------------------------------------------

import scanpy as sc # single-cell analysis toolkit
import os
import pandas as pd  # For DataFrame creation and CSV/TSV writing

#from run_scanpy_paga_fdg_clusters import sample_id

#-----------------------------------------------------------------------------------
#   Params
#-----------------------------------------------------------------------------------

# set PAGA group
# paga_group = "celltypes"
# paga_group = "clusters"

# set the sample ID
#sample_id = "e11_control"
#sample_id = "e11_ko"
#sample_id = "e12_control"
sample_id = "e12_ko"

# set the main project directory
project_dir = f"/home/melanie-smith/workDir/sophieWiszniak/20251104_sophieWiszniak_ncc_pa"

# set the loom file address
loom_file = os.path.join(
    project_dir,
    "outDir",
    "09-rna_velocity",
    sample_id,
    "velocyto_output",
#    "possorted_genome_bam_5Q1UF.loom" # e11_control
#    "possorted_genome_bam_YIL7N.loom" # e12_control
#    "possorted_genome_bam_LZRTJ.loom" # e11_ko
    "possorted_genome_bam_0YX73.loom" # e12_ko
)

# set the output directory
out_dir = os.path.join(
    project_dir,
    "outDir",
    "09-rna_velocity",
    sample_id,
    "velocyto_output"
    )

os.makedirs(out_dir, exist_ok=True)

#-----------------------------------------------------------------------------------
# Import the .loom file (spliced/unspliced counts)
#-----------------------------------------------------------------------------------

# import the data from the 10x velocyto run (counts are spliced, unspliced and ambiguous)
vlm = sc.read_loom(loom_file)


# Extract the spliced layer as a sparse matrix.
spliced_matrix = vlm.layers['spliced']

# Convert the sparse matrix to a dense NumPy array.
# Note: This step densifies the matrix, which for your dimensions (9061 cells x 32286 genes)
# will require approximately 2.3 GB of RAM. If memory is an issue, consider exporting
# in a sparse format (e.g., Matrix Market .mtx) instead, but since you requested .tsv/.csv,
# we'll proceed with dense.
spliced_dense = spliced_matrix.toarray()

# Create a Pandas DataFrame from the dense array.
# Use cell names (obs_names) as row index and gene names (var_names) as column headers.
# This structure matches typical count matrix formats for import into Seurat.
spliced_df = pd.DataFrame(
    spliced_dense,
    index=vlm.obs_names,
    columns=vlm.var_names
)

# Write the spliced DataFrame to a compressed TSV file.
# We use .tsv.gz for efficiency (gzip compression reduces file size significantly,
# as the matrix is sparse with many zeros). In R, you can read this with data.table::fread()
# or readr::read_tsv(), which handle .gz files natively.
# float_format='%.0f' ensures integer-like output if values are whole numbers (common for counts).
spliced_df.to_csv(
    os.path.join(out_dir, "spliced_counts.tsv.gz"),
    sep='\t',
    index=True,      # Include cell names as the first column
    header=True,     # Include gene names as the header row
    compression='gzip',
    float_format='%.0f'
)

# Step 6: Repeat the process for the unspliced layer.
unspliced_matrix = vlm.layers["unspliced"]
unspliced_dense = unspliced_matrix.toarray()
unspliced_df = pd.DataFrame(
    unspliced_dense,
    index=vlm.obs_names,
    columns=vlm.var_names
)
unspliced_df.to_csv(
    os.path.join(out_dir, "unspliced_counts.tsv.gz"),
    sep='\t',
    index=True,
    header=True,
    compression='gzip',
    float_format='%.0f'
)

# Print confirmation messages
print("Exported spliced counts to 'spliced_counts.tsv.gz'")
print("Exported unspliced counts to 'unspliced_counts.tsv.gz'")

