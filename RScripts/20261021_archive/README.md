# Harmony batch-integration check — archived 2026-10-21

## Why this exists

Sophie asked about transitions within the Mesenchyme population. Before
looking at that, we revisited an earlier, never-quantified assumption that
batch effects weren't an issue in this dataset (the original check was just
a visual PCA/UMAP coloured by `orig.ident` in `02-initial_QC.Rmd` — no
mixing metric was ever computed). This archive is a self-contained side
exploration of whether Harmony integration changes that picture, and what
that means for any trajectory/transition work on Mesenchyme. It was never
folded into the main numbered pipeline and doesn't replace
`final_interesting_cells_0.6.rds` or any other active pipeline object.

**Status: parked, not resolved.** Melanie is discussing options with Sophie
before deciding whether/how to bring this back into the main pipeline.

## What's in here

| File | What it does |
|---|---|
| `06.5-harmony_batch_integration_check.Rmd` | Loads `final_interesting_cells_0.6.rds`, splits the RNA layer by `stage` (batch = sequencing run, e11 vs e12 — confirmed as a deliberate choice over splitting by individual sample), runs `IntegrateLayers(method = HarmonyIntegration)`, clusters/embeds at the same parameters as the original (`dims = 1:30`, `resolution = 0.6`), and produces original-vs-Harmony comparison plots. |
| `06.5b-mesenchyme_substructure_harmony_fork.Rmd` | Joins the mesenchyme subtype labels (from `PScripts/14.9-export_subtype_metadata_for_propeller.py`'s output CSV) onto the Harmony-integrated object and plots them on both the original and Harmony UMAPs, to see how the EndoMT/E11_NCC/E12_NCC sub-structure looks before vs after integration. |
| `06.5c-export_embeddings_for_fdg_fork.Rmd` | Exports the PCA and Harmony cell embeddings + metadata to CSV, for the PAGA/FDG check below. |
| `exploratory_mesenchyme_vsmc_fdg_fork.py` | (Originally in `PScripts/` — moved here to keep this whole piece of work together.) Runs PAGA + a force-directed graph on both embeddings, same recipe/parameters as `PScripts/14.3-run_scanpy_paga_fdg_cell.py` (pcs=30, neighbours=25), to check whether an apparent UMAP gap between Mesenchyme and VSMC reflects a real discontinuity. |

None of these scripts modify any existing pipeline script or object — each
only reads already-saved outputs (`final_interesting_cells_0.6.rds`, the
`14.9` subtype CSV) and writes to its own output folder.

## Where the output went

Everything these scripts produced is archived alongside, in
`outDir/20261021_archive/06.5-harmony_batch_integration_check/` (moved from
its original location at `outDir/06.5-harmony_batch_integration_check/`):

- `final_interesting_cells_0.6_harmony_integrated.rds` — the Harmony-integrated
  Seurat object (adds `harmony`/`umap.harmony` reductions and
  `harmony_clusters`, alongside the untouched original `pca`/`umap`/
  `seurat_clusters`)
- `compare_umap_by_stage.png`, `compare_umap_by_genotype.png`,
  `compare_umap_by_celltype.png` — original vs Harmony, full dataset
- `compare_umap_mesenchyme_subtype.png`,
  `compare_umap_mesenchyme_subtype_by_stage.png` — mesenchyme subtype labels
  on both embeddings
- `cluster_crosstab_original_vs_harmony.csv`,
  `celltype_by_harmony_cluster.csv` — cluster reassignment cross-tabs
- `mesenchyme_vsmc_fdg_check/` — PAGA connectivity matrices, FDG plots, and
  the underlying `.h5ad`/CSV exports for the Mesenchyme-VSMC bridge check

## What we found (summary for picking this back up later)

1. **The earlier "no batch effect" call doesn't hold up for Mesenchyme.**
   On the original (uncorrected) UMAP, a large Mesenchyme sub-population
   separates almost entirely by stage (e11 vs e12). Harmony integration
   resolves that separation.

2. **Genotype is not confounded with batch** — well-mixed before and after
   integration, in both the full-dataset view and the subtype-level view.
   Genotype-driven findings elsewhere in the project aren't called into
   question by any of this.

3. **Stage and batch are perfectly confounded in this dataset** (every e11
   sample was sequenced separately from every e12 sample), so Harmony
   cannot distinguish "technical noise" from "real e11→e12 developmental
   change." Whatever it aligns along that axis, it aligns indiscriminately.

4. **The stage-driven separation is not uniform across cell types** — this
   is the key finding for deciding what to do next. Quantifying "% of a
   cell's nearest neighbours that are the same stage" vs the random-mixing
   expectation for each cell type:

   | Cell type | Observed same-stage neighbours | Expected under random mixing | Excess |
   |---|---|---|---|
   | Mesenchyme | 87.5% | ~50.5% | **+37.0** |
   | VSMC | 87.9% | ~59.7% | **+28.2** |
   | Myocardium | 75.4% | ~55.1% | +20.3 |
   | Epicardium | 75.9% | ~55.1% | +20.8 |
   | Endocardium | 72.8% | ~50.6% | +22.2 |

   Every cell type shows some excess stage segregation (~+20 points),
   consistent with a real background technical effect touching the whole
   dataset roughly evenly — which is exactly what you'd expect from a
   genuine sequencing-run batch effect, since all cells from the same
   library go through the same prep. But Mesenchyme (+37) and VSMC (+28)
   sit well above that common baseline, with Mesenchyme the biggest outlier
   by a clear margin.

   A pure technical batch effect has no obvious reason to hit two specific
   cell types harder than the other three. But Mesenchyme and VSMC are
   exactly the two neural-crest-derived lineages known to be actively
   differentiating during this exact developmental window: NCC-derived
   mesenchyme is transitioning/maturing, and VSMC differentiation of cardiac
   NCC derivatives in the outflow tract/great vessels is one of the
   best-characterized NCC fate decisions happening right around
   E11.5-E12.5. Endocardium, Epicardium, and Myocardium, by contrast, are
   comparatively more terminally specified/structurally stable cell types by
   this stage, so there's less genuine biological turnover between e11 and
   e12 for them — leaving mostly the shared technical baseline.

   **Take-away:** the data support a layered explanation rather than a
   either/or one — a uniform ~20-point technical effect present everywhere,
   plus a real, additional biological signal concentrated in the two
   lineages that are genuinely changing state fastest during this window.
   This is the opposite of reassuring for Harmony-correcting Mesenchyme
   specifically: it means a meaningful share of what gets "corrected away"
   there is likely real developmental signal, not just batch noise — so
   treating the Harmony-integrated embedding as a straightforward fix for
   Mesenchyme is not safe to do without flagging this explicitly.

5. **The apparent Mesenchyme/VSMC gap in UMAP is a rendering artifact, not a
   real discontinuity** — unrelated to the batch question above. PAGA
   connectivity and a force-directed graph (same neighbour graph, no UMAP
   involved) show Mesenchyme and VSMC directly adjacent/touching in both the
   original and Harmony embeddings, with a connectivity score 3-4x stronger
   than the already-characterized (and already-trusted) 436-cell
   Endocardium-Mesenchyme bridge.

## Open question for Sophie

Given stage and batch are perfectly confounded, is it better to (a) use the
Harmony-integrated embedding for mesenchyme transition work, accepting that
some real biology may be smoothed along with batch noise, (b) stay with the
original embedding and treat the stage-driven separation as a caveat rather
than correct it, or (c) a middle path — regress out only known technical
covariates (depth, % mito) without aligning on stage directly, keeping stage
as an explicit modelled covariate downstream instead. No option is fully
clean with this experimental design (stage and batch cannot be separated
after the fact); see the full discussion in the chat session from
2026-10-21 for the detailed pros/cons.
