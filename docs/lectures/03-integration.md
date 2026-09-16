# Lecture 3 — Integration

ngs101 Part 3. Scope: per-sample QC'd matrices in, **one merged, batch-corrected
object** out. Clustering deliberately belongs to the next stage.

Setup: 8 samples, 4 healthy and 4 post-treatment. The tutorial tries CCA, RPCA,
Harmony and FastMNN, and scores them with silhouette and local mixing.

---

## 1. The goal, and how to know you have a problem

Merge 8 samples, run PCA and UMAP, colour by sample. If cells group by sample,
something is wrong — you want cells organised by **cell type**, with condition
appearing as variation *within* those types, and donor explaining almost nothing.

But sample-driven separation is not automatically a batch effect. Treatment can
genuinely remove a population or shift a state, and "correcting" that deletes the
result. Four ways to tell them apart:

1. **Are the same cell types split by sample?** If CD14 monocytes from two donors
   form separate islands but both express `CD14`/`LYZ`, it is one cell type
   shifted — batch. A group present only in some samples, with its own markers,
   may be real.
2. **Does separation follow processing or condition?** Splitting by collection day
   regardless of condition is technical. Treated consistently separating from
   healthy across different donors *and* days is plausibly biology.
3. **Global or cell-type-specific?** Technical effects — depth, ambient, stress —
   push *every* cell type the same direction. Biology is usually selective
   (interferon response in treated monocytes but not B cells).
4. **Which genes load the PCs?** `FOS`/`JUN`/`HSPA1A`, mito, ribo, or
   depth-correlated genes → technical. `ISG15`/`IFI6` only in treated → biology.

## 2. Confounding — the failure no algorithm fixes

If all 4 healthy samples were processed on day 1 and all 4 treated on day 2, then
**batch and condition are the same variable**. Any correction that removes the
day effect removes the treatment effect, and nothing in the data can separate
them. This is an experimental-design failure; the only fix is before sequencing.

`integrate.py` refuses to run when `batch_key` is perfectly confounded with
`condition_key`, rather than emitting a corrected object with the biology
quietly removed.

## 3. Over-correction

Every method assumes batches share most cell types. When one condition has a
population the others lack, an aggressive method forces those cells onto the
nearest shared cluster and a treatment-specific population disappears into its
neighbour.

So the choice of method is never blind: you score **mixing** (are samples
interleaved?) against **conservation** (are cell types still distinct?). Both
matter, and they pull against each other.

## 4. You must declare the batch variable

Nothing infers it.

| Tool | How |
|---|---|
| Harmony (R) | `RunHarmony(obj, group.by.vars = "sample_id")` |
| Harmony (Py) | `sc.external.pp.harmony_integrate(adata, key = "sample")` |
| Seurat v5 | split layers by a metadata column, then `IntegrateLayers()` |
| FastMNN | `batch =` |
| scVI | `batch_key = "sample"` |

- **Per sample** is the safe default when each was prepped and loaded separately.
- **Per run or chip** is more correct when several samples were processed
  together — correcting per-sample there removes genuine donor biology.

Harmony accepts several (`group.by.vars = c("donor","chemistry")`), but each one
removes more variation; over-specifying is a real way to over-correct.

**Never pass condition as the batch key.** That removes what you are studying.

This is what the samplesheet's `batch` column has been for since Part 1:
`meta.batch = row.batch ?: row.sample`.

## 5. Anchors, and what has to be true for them to exist

Harmony and FastMNN work on the PCA embedding. CCA and RPCA are Seurat's
**anchor-based** methods.

**An anchor is a *pair of cells*, one from each sample** — not genes, not a
matrix. Take cell *a* in sample 1 and cell *b* in sample 2: if *b* is among *a*'s
nearest neighbours **and** *a* is among *b*'s (mutual nearest neighbours), they
are probably the same cell type in two samples. Each pair yields a difference
vector — "this is how far sample 2 is shifted for this kind of cell" — and their
average is the correction applied to everything else.

Shared genes define the space the anchors are found in (`SelectIntegrationFeatures`);
anchors are cells within it.

- **CCA** finds gene combinations maximally correlated between datasets. Powerful
  when the batch effect is large, but it *insists* on finding shared structure,
  so it over-corrects more readily.
- **RPCA** projects each dataset into the other's PCA space. Conservative, much
  faster, right when samples share less biology or datasets are large. Seurat's
  own guidance: reach for RPCA when CCA looks like it is over-merging.

Seurat also *scores* anchors — checking whether the two cells' neighbourhoods
overlap — and discards low-scoring ones, a filter against spurious pairs.

**The assumption: samples must share cell populations.** Anchors only exist if
the same types are present in both. MNN will still pair a unique population with
*something*, because the algorithm cannot say "no match exists" — which is
exactly the over-correction failure above.

FastMNN is the same MNN idea applied directly in PCA space.

## 6. What integration actually modifies

**An embedding, not your expression data.** Harmony takes the PCA coordinates
(~30 numbers per cell) and returns corrected coordinates. Counts and
log-normalised values are never touched.

Harmony's loop, until convergence: soft-cluster in PC space with a penalty
rewarding batch-mixed clusters → compute each batch's centroid within each
cluster → shift cells toward the shared centre, weighted by cluster membership.

| Method | Returns |
|---|---|
| Harmony | corrected PCA embedding |
| Seurat v5 `IntegrateLayers` | corrected reduction, e.g. `integrated.cca` |
| FastMNN | corrected low-dimensional reduction |
| scVI | latent space |

**Seurat v4 was different** — `IntegrateData()` produced an `integrated` assay of
corrected *expression* values, and much published work wrongly ran DE on it. v5
returns reductions partly to stop that.

### Why this matters for DE

1. **Integration deliberately removes between-sample variation** — which DE needs
   to estimate uncertainty. Remove it and p-values become falsely significant.
2. **It can remove the condition effect**, which lives between samples.
3. **Corrected values are not counts.** DESeq2/edgeR model counts with a negative
   binomial; shifted continuous values break the model.

Donors belong in the DE *model* — typically pseudobulk per donor with condition
as the tested variable (Part 5).

### Our output contract

```
layers["counts"]      raw counts, untouched     -> DE (Part 5)
X                     log-normalised, untouched -> marker plots
obsm["X_pca"]         uncorrected PCA
obsm["X_integrated"]  corrected embedding       -> neighbours, clustering, UMAP
```

Nothing downstream ever reads "corrected expression".

## 7. Scoring: two axes, and why local mixing earns its place

**Silhouette width** per cell: `s = (b − a) / max(a, b)`, where `a` is mean
distance to its own group and `b` to the nearest other group.

- by **cell type** → high is good (biology survived)
- by **sample** → low is good (batch removed)

**Local mixing**: for each cell, how diverse are the sample labels among its 50
nearest neighbours (entropy or inverse Simpson — the idea behind **iLISI**).

Two silhouettes *do* cover both axes — local mixing is not a third axis but a
different **resolution** on the batch axis. It earns inclusion for three reasons:

1. **Silhouette is global, mixing is local.** The commonest residual batch effect
   is per-sample sub-islands *within* each cell type: samples overlap globally
   while every cell's 50 neighbours remain same-sample. Silhouette is the metric
   least able to see this.
2. **Silhouette assumes groups are blobs.** "Batch" is not — sample 3's cells are
   scattered across every cell type, so mean distance to a multimodal group is a
   shaky statistic. kNN mixing assumes nothing.
3. **It matches what actually breaks.** Clustering runs on the kNN graph. If
   neighbours are same-sample, Leiden carves out sample-specific clusters
   regardless of the global silhouette.

Practical: silhouette needs all pairwise distances, so the tutorial subsamples to
5,000 cells, weakening it further; kNN mixing runs on every cell.

**Each metric alone has a degenerate optimum.** Maximise mixing by collapsing
everything into one blob (perfect over-correction scores perfectly); maximise
conservation by doing nothing (perfect under-correction scores perfectly). This
is why scIB reports both and why the honest result is a trade-off, not a winner.

**Our caveat:** true cell-type labels do not exist yet — annotation is Part 4. We
use clusters from the **unintegrated** embedding as a proxy, which is itself
partly batch-driven. It catches gross over-correction but is not ground truth,
and `integration_qc.py` says so in its output.

## 8. Default choice

**Harmony**, because:

- it only touches the embedding (~30 numbers/cell), so it scales
- counts stay pristine for Part 5
- it benchmarks near the top in scIB for batch removal with decent conservation
- the *same algorithm* exists in both languages (`harmony` / `harmonypy`), so it
  satisfies the "cross-language only when genuinely good" rule without
  implementing two different methods
- multiple covariates come free

**Switch when:**

| To | When |
|---|---|
| RPCA | Harmony looks over-corrected — an expected condition-specific population vanished |
| CCA | batch effect huge relative to biology (cross-technology, cross-species) |
| scVI/scANVI | large atlases, probabilistic latent space, label transfer |
| none | a single batch |

That last one is a guard, not a setting: with one `batch_key` level the pipeline
**skips integration and logs it**, rather than running Harmony on one batch and
producing something meaningless.

Limitation to hold onto: because Harmony returns only an embedding, it cannot
produce corrected expression at all. For us that is a feature — and if a
downstream tool demands corrected expression, the honest answer is usually that
it should not.
