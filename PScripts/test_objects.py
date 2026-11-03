import scvelo as scv
import scanpy as sc
import os
import numpy as np


#-----------------------------------------------------------------------------------
#   Params
#-----------------------------------------------------------------------------------

# set PAGA group
paga_group = "celltypes"


# set the sample ID
sample_id = "e11_control"

# set the number of PCs
pcs="pcs20"


# set the number of neighbours
neighbours="neighbours15"


# set the count level (ie gene or spliced/unspliced)
count_level = "spliced_unspliced"

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
    f"{sample_id}_fdg_{pcs}_{neighbours}_counts_{count_level}_{paga_group}.h5ad")

# ── 1. Load both files ────────────────────────────────────────────────────────
adata = sc.read_h5ad(h5ad_file_input)
vlm   = sc.read_loom(loom_file)

print("=" * 60)
print("ADATA (h5ad)")
print("=" * 60)
print(f"Shape:          {adata.shape}")
print(f"obs_names[0]:   {adata.obs_names[0]}")
print(f"Layers:         {list(adata.layers.keys())}")
print(f"X dtype:        {adata.X.dtype}")
print(f"obsm keys:      {list(adata.obsm.keys())}")
print(f"uns keys:       {list(adata.uns.keys())}")
print()

print("=" * 60)
print("VLM (loom)")
print("=" * 60)
print(f"Shape:          {vlm.shape}")
print(f"obs_names[0]:   {vlm.obs_names[0]}")
print(f"Layers:         {list(vlm.layers.keys())}")
print()

# ── 2. Normalise barcodes (same as your script) ───────────────────────────────
adata.obs_names = adata.obs_names.str.split("_").str[0]
vlm.obs_names   = vlm.obs_names.str.split(":").str[1].str.rstrip("x") + "-1"

# ── 3. Overlap diagnostics ────────────────────────────────────────────────────
common     = adata.obs_names.intersection(vlm.obs_names)
only_adata = adata.obs_names.difference(vlm.obs_names)
only_vlm   = vlm.obs_names.difference(adata.obs_names)

print("=" * 60)
print("BARCODE OVERLAP")
print("=" * 60)
print(f"adata cells:        {len(adata)}")
print(f"loom cells:         {len(vlm)}")
print(f"Common:             {len(common)}")
print(f"Only in adata:      {len(only_adata)}")
print(f"Only in loom:       {len(only_vlm)}")

if len(only_adata) > 0:
    print(f"\n  adata-only examples: {only_adata[:5].tolist()}")
if len(only_vlm) > 0:
    print(f"  loom-only  examples: {only_vlm[:5].tolist()}")

# ── 4. Check spliced/unspliced layers exist in loom ───────────────────────────
print()
print("=" * 60)
print("LAYER CONTENT CHECK")
print("=" * 60)
for layer in ["spliced", "unspliced", "ambiguous"]:
    if layer in vlm.layers:
        mat = vlm.layers[layer]
        total = mat.sum()
        print(f"  loom['{layer}']:  shape={mat.shape}, total counts={total:.0f}")
    else:
        print(f"  loom['{layer}']:  *** MISSING ***")

# ── 5. Check adata for conflicting spliced layer ──────────────────────────────
print()
print("=" * 60)
print("ADATA LAYER CONTENT CHECK")
print("=" * 60)
for layer in list(adata.layers.keys()):
    mat = adata.layers[layer]
    print(f"  adata['{layer}']: shape={mat.shape}, total={mat.sum():.0f}")

# ── 6. Attempt merge and inspect result ──────────────────────────────────────
print()
print("=" * 60)
print("POST-MERGE CHECK")
print("=" * 60)
vlm_sub = vlm[common].copy()
vlm_sub = vlm_sub[adata.obs_names].copy()
merged  = scv.utils.merge(adata, vlm_sub)

print(f"Merged shape:   {merged.shape}")
print(f"Merged layers:  {list(merged.layers.keys())}")
for layer in ["spliced", "unspliced"]:
    if layer in merged.layers:
        mat = merged.layers[layer]
        print(f"  merged['{layer}']: shape={mat.shape}, nnz={mat.nnz if hasattr(mat,'nnz') else 'dense'}, total={mat.sum():.0f}")
    else:
        print(f"  merged['{layer}']: *** MISSING — this will break scv.pl.proportions ***")

# ── Gene name format check ────────────────────────────────────────────────────
print("ADATA var_names examples:")
print(adata.var_names[:10].tolist())

print("\nLOOM var_names examples:")
print(vlm.var_names[:10].tolist())

# Check for Ensembl IDs vs gene symbols
adata_ensembl = adata.var_names.str.startswith("ENSM").sum()
vlm_ensembl   = vlm.var_names.str.startswith("ENSM").sum()
print(f"\nadata Ensembl IDs: {adata_ensembl} / {adata.n_vars}")
print(f"loom  Ensembl IDs: {vlm_ensembl}  / {vlm.n_vars}")

# Check overlap
shared = adata.var_names.intersection(vlm.var_names)
print(f"\nShared genes (exact match): {len(shared)}")

# If overlap is low, check if loom has both gene_id and gene_name columns
print("\nLoom var columns:")
print(vlm.var.columns.tolist())
print(vlm.var.head(10))

print("\nAdata var columns:")
print(adata.var.columns.tolist())
print(adata.var.head(10))