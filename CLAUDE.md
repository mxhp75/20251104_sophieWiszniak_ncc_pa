# CLAUDE.md — Project Instructions for Claude Code

## Behaviour Rules — Read Before Acting

- **Always read a file in full before editing it**
- **Do not make changes to any file without first describing what you plan to change and waiting for confirmation**
- When in doubt about intent, ask — do not guess and proceed
- Prefer minimal, targeted edits over rewrites
- Do not add new functionality, parameters, or restructure code unless explicitly asked

---

## Project Overview

**Project:** Neural Crest Cell (NCC) Pharyngeal Arch scRNA-seq analysis: Mib1 knock out
**Researcher:** Melanie Smith — melanie.smith@adelaide.edu.au  
**PI/Collaborator:** Sophie Wiszniak  
**Project directory:** `/home/melanie-smith/workDir/sophieWiszniak/20251104_sophieWiszniak_ncc_pa`

This project analyses single-cell RNA sequencing data from mouse embryos (E11.5 and E12.5) with a knockout (KO, Mib1) and control condition. The biological focus is on neural crest cell lineage trajectories in the developing pharyngeal arch, including differentiation into Myocardium, VSMC, Endocardium, Epicardium, and Mesenchyme cell types. Variance partition was run prior to this workflow to establish genes who's variance was associated with cell cycle. The genes with the highest variance were removed  prior to the current workflow.

**Known technical caveat — `Mib1` read counts cannot validate the knockout**: the model is `Wnt1-Cre; Mib1^fl/fl`, a conditional knockout where Cre excises an *internal* exon of `Mib1` from genomic DNA, only in neural-crest-lineage cells. The gene's promoter, 5′ end, and (typically) final exon/3′UTR/poly-A site are untouched. 10x 3′ scRNA-seq only sequences the 3′ end of each transcript, so a recombined (non-functional) transcript looks identical to normal `Mib1` from that vantage point — this assay cannot detect the deletion. Confirmed directly (2026-10-20, against `final_interesting_cells_0.6.rds`): `Mib1` expression is not reduced in KO — if anything slightly higher — consistently across every sample, stage, and cell type. This is expected given the above, not evidence against a successful knockout. Don't use `Mib1` read counts in this data to assess knockdown efficiency; the EGFP reporter (if linked to Cre recombination) is a better proxy for which cells recombined. Exact exon(s) targeted in this `Mib1^fl` allele and the Cre/EGFP reporter design are unconfirmed with Sophie as of this date. Full writeup: `outDir/14-rna_velocity/additional_plots_meeting_20261013_lab_book.md`, section 7.

---

## Samples

Four samples, each processed independently:

| Sample ID     | Stage | Condition |
|---------------|-------|-----------|
| `e11_control` | E11   | Control   |
| `e11_ko`      | E11   | Knockout  |
| `e12_control` | E12   | Control   |
| `e12_ko`      | E12   | Knockout  |

---

## Workflow — Step 10: RNA Velocity

All R and Python scripts are in `/home/melanie-smith/workDir/sophieWiszniak/20251104_sophieWiszniak_ncc_pa/Scripts`, numbered as one shared pipeline sequence (language is indicated only by file extension). Outputs write to:
`outDir/10-rna_velocity/{sample_id}/`

### Scripts in order

| Script | Language | Description |
|--------|----------|-------------|
| `10.1-subset_seurat_orig_ident.Rmd` | R | Splits the full merged Seurat object into one object per sample and saves as `{sample_id}_no_cc.rds` |
| `10.2-process_seurat_samples.Rmd` | R | Re-embeds each sample-level Seurat object (PCA, UMAP, clustering) and saves as `{sample_id}_no_cc_re_embed.rds` |
| `10.3-loom_to_seurat_samples.Rmd` | R | Imports Velocyto spliced/unspliced counts from `.loom` files, builds a combined Seurat object, and exports as `{sample_id}_spliced_unspliced.h5ad` |
| `10.1-run_scanpy_paga_fdg_cell.py` | Python | Reads the `.h5ad` files, runs PAGA and Force Directed Graph (FDG) across a grid of PC and neighbour parameters, saves plots and `.h5ad` files to `scanpy_output/celltypes/` |
| `10.2-run_rnaVelocity_dynamical_fdg_celltypes.py` | Python | Runs scVelo dynamical RNA velocity on the PAGA/FDG `.h5ad` files, merging in the Velocyto `.loom` data. Loops over all four samples. |

### Key file naming conventions

- PAGA/FDG `.h5ad` input to velocity script:
  `outDir/10-rna_velocity/{sample_id}/scanpy_output/celltypes/{sample_id}_fdg_pcs{n_pcs}_neighbours{n_neighbours}_counts_{count_level}_celltypes.h5ad`

- Velocyto `.loom` files live in:
  `outDir/09-rna_velocity/{sample_id}/velocyto_output/`

- Loom barcode prefix per sample:

| Sample        | Loom prefix                        |
|---------------|------------------------------------|
| `e11_control` | `possorted_genome_bam_5Q1UF:`      |
| `e11_ko`      | `possorted_genome_bam_LZRTJ:`      |
| `e12_control` | `possorted_genome_bam_YIL7N:`      |
| `e12_ko`      | `possorted_genome_bam_0YX73:`      |

### Barcode format notes

- Velocyto loom barcodes: `{prefix}:{barcode}x` (e.g. `possorted_genome_bam_5Q1UF:AAATGGATCGATACGTx`)
- Seurat/scanpy adata barcodes: `{barcode}-1_{suffix}` (e.g. `AAATGGATCGATACGT-1_1_1`)
- Cleaning steps: strip prefix on `:`, replace trailing `x` with `-1`, append Seurat suffix (auto-detected from `adata.obs_names`)

---

## Environment

- R environment: standard Bioconductor/Seurat stack
- Python environment: `scv` conda env (`miniconda3/envs/scv`), Python 3.9
- Key Python packages: `scvelo`, `scanpy`, `matplotlib`
- Note: `scv.tl.recover_dynamics()` must be run with `n_jobs=1` in interactive/IPython sessions due to a `BrokenPipeError` with the `loky` multiprocessing backend. Use `n_jobs=4` when running from the terminal.
