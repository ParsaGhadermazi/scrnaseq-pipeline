# Lecture 1 -- FASTQ to count matrix

Scope: raw reads in, count matrices out. Part 2 starts from what this produces.

---

## 1. Where per-cell identity comes from

Bulk RNA-seq gives one column per library. 10x gives ~10,000 columns from a
**single** pooled library -- there was never a tube per cell. The identity is
written into the molecules themselves, before pooling.

**The mechanism.** The Chromium controller emulsifies cells into droplets. Each
droplet -- a **GEM** (Gel bead-in-EMulsion) -- receives one gel bead coated with
millions of copies of one oligo. All those copies share a single **16 bp cell
barcode**, but each carries a different **12 bp UMI**.

Inside the droplet: the cell lyses, the bead dissolves, and reverse transcription
runs. The barcode and UMI become covalently part of the cDNA. Only then is the
emulsion broken and everything pooled -- safe, because identity is now sequence.

The sequencer has no idea it is sequencing single cells. Demultiplexing is
entirely computational, done afterwards by reading barcodes off R1.

## 2. Barcode collisions

The v3 whitelist holds ~6.8M distinct barcode *sequences*, but 10x manufactures
*billions* of beads. Many beads share a sequence, and beads are drawn at random.
That is the birthday problem:

    expected collisions = k^2 / 2N = 10,000^2 / (2 * 6.8e6) ~ 7 pairs

~7 of 10,000 cells silently merge with another. The alternative -- a unique
barcode per bead -- would put barcodes one base apart in sequence space, so a
single miscalled base would reassign a read to the wrong cell. **10x trades ~7
collisions for unambiguous error correction.** A curated whitelist also lets
Cell Ranger reject invalid barcodes and repair 1-base errors.

A collided barcode is two transcriptomes summed -- indistinguishable from a
doublet, and part of what Part 2's doublet detection removes.

## 3. The UMI: molecules, not reads

A cell holds ~1e5-5e5 mRNA molecules -- too little to sequence, so it is
PCR-amplified. **PCR is exponential and uneven**: one transcript amplifies 500x,
its neighbour 50x, for reasons unrelated to abundance. Counting reads would
largely measure amplification luck.

So each molecule is tagged with a random 12 bp UMI *before* amplification. Every
copy inherits it. Then:

    count of DISTINCT UMIs for (cell barcode, gene) = original molecule count

4^12 ~ 16.7M possible tags, so same-gene same-cell collisions are rare -- though
at high expression they do happen (**UMI collision**), slightly undercounting the
most-expressed genes.

Related metric: **sequencing saturation**. Sequence deeper and you stop finding
new UMIs and merely re-read known ones. High saturation means more sequencing
buys nothing.

## 4. Droplet loading is Poisson

You cannot aim a cell into a droplet. With lambda = cells / droplets:

| cells in droplet | P (lambda = 0.1) | meaning |
|---|---|---|
| 0 | 90.5% | empty droplet |
| 1 | 9.0%  | what you want |
| >=2 | 0.47% | multiplet |

~100,000 droplets are formed to recover ~10,000 cells. **90% are deliberately
wasted** -- emptiness is the price of keeping multiplets rare.

As a fraction of *occupied* droplets, multiplets are 0.47/9.47 ~ **4.9%**, and for
small lambda that is almost exactly **lambda/2**. So the doublet rate scales
*linearly* with cells recovered -- which is where 10x's **0.8% per 1,000 cells
recovered** rule comes from. 5,000 cells -> ~4%; 10,000 -> ~8%.

lambda is not directly measurable: capture efficiency is only ~50-65%, so
recovering 10,000 cells means loading ~16,000-20,000.

**The trade.** Load denser -> more cells, proportionally more doublets. Load
sparser -> cleaner, but a channel costs the same either way, and a 1%-frequency
cell type gives 100 cells at 10,000 recovered versus 20 at 2,000, often too few
to cluster at all. Doublets are not avoidable, they are **priced**.

The chip's ~80-100k GEM capacity is fixed, so the only real knob is cell
concentration.

## 5. How doublet detection uses that number

You cannot train on known doublets -- there are no labels. So you manufacture them:

1. Sum two random real cells -> a synthetic doublet, known fake (`pN = 0.25`
   means synthetics are 25% of the pool).
2. Embed real + synthetic together in PCA space.
3. For each *real* cell: what fraction of its *k* nearest neighbours are
   synthetic? That fraction is **pANN**. `k` is the parameter **pK**.

**DoubletFinder** then requires `nExp` -- the expected doublet *count*, computed
by hand from the 0.8%-per-1,000 rule -- and calls exactly the top `nExp` cells.
It also needs `paramSweep()` -> `find.pK()` to tune `k` per dataset.

**scDblFinder** estimates the rate itself and picks its own threshold: no `nExp`,
no sweep. That is why Part 2-2 makes it the default. **Scrublet** and
**DoubletDetection** are the Python implementations of the same idea.

## 6. Cell calling: knee plots and EmptyDrops

Empty droplets are not empty of RNA -- the suspension contains free mRNA from
cells lysed during dissociation (**the ambient soup**), so every barcode has some
counts.

**Knee plot.** Sort barcodes by UMI count, plot rank vs count log-log: a high
plateau (cells), a sharp cliff (the knee), a long tail (ambient). It is the first
figure in `web_summary.html`, and a soft, smeared knee means the dissociation
went badly. Early Cell Ranger thresholded here directly.

**Its failure mode:** genuinely small cells -- platelets, neutrophils, quiescent
lymphocytes -- carry little RNA, land in the ambient range, and get discarded.
You systematically lose the hardest-to-study cell types.

**EmptyDrops** fixes this by testing *profile*, not depth:

- Pool barcodes below a low cutoff (`lower = 100`) into an **ambient profile**.
- For each candidate: is its profile significantly different from ambient, or
  just a bigger scoop of the same soup? Dirichlet-multinomial likelihood,
  Monte Carlo p-value (`niters = 10000`).
- Keep `FDR < 0.01`.

A low-count barcode now survives if its transcriptome is *distinctive*. Cell
Ranger v3+ combines both: obvious high-count barcodes first, EmptyDrops below.

**The ambient profile EmptyDrops estimates is the same soup SoupX subtracts in
Part 2.** Both need the empty droplets -- which is why the pipeline must carry
the **raw** matrix forward, not just the filtered one.

Python counterpart: **CellBender**, which does cell calling and ambient removal
jointly with a generative model.

## 7. The reference package

`cellranger count` refuses a bare genome FASTA. It needs a prebuilt package
(~11 GB): genome, STAR index, and a **filtered** GTF.

**Why filtered.** 10x keeps only certain biotypes (protein_coding, lncRNA, IG/TR
segments) and strips **pseudogenes**, readthrough transcripts, and retained-intron
variants. The rule that makes this matter:

> Cell Ranger counts a read only if it maps uniquely to **one** gene. A read
> overlapping two annotated features is discarded as ambiguous.

So a pseudogene resembling its parent doesn't merely add a junk row -- it
**silently deletes counts from the real gene**.

**The intron decision.** ~20-30% of 10x reads are intronic: the poly-dT primer
binds A-rich stretches inside introns, and droplets capture nuclear pre-mRNA.
For single-*nucleus* work, most RNA is unspliced. 10x switched the default to
`--include-introns=true` in **Cell Ranger 7.0**. Counts rise substantially, which
shifts UMI distributions, which shifts cell calling, which shifts everything.

**Two runs differing only in this flag are not comparable** -- so the pipeline
exposes it *and* records it.

Consequence for later: spliced vs unspliced separation is the entire basis of
**RNA velocity**. Cell Ranger merges them into one number and the distinction is
gone, which is why Part 13 (scVelo) needs `velocyto` or an alevin-fry splici
index. A decision at step one determines whether Part 13 is possible at all.

## 8. Tooling: SRA, file layout, naming

**Read structure (3' v3):**

```
P5-[i5]-[Read1 primer]-[16bp CELL BARCODE][12bp UMI][polydT]-[ cDNA ]-[Read2 primer]-[i7]-P7
                        |________ R1 = 28 bp ________|        |R2=91bp|
```

| File | Length | Contents | Identifies |
|---|---|---|---|
| R1 | 28 bp | 16 bp barcode + 12 bp UMI | which **cell**, which **molecule** |
| R2 | 91 bp | cDNA near the transcript 3' end | which **gene** |
| I1 | 8 bp  | i7 sample index | which **library** |
| I2 | 10 bp | i5 sample index (dual-index) | which **library** |

**Two separate levels of demultiplexing:** sample-level via I1/I2, done by
`bcl2fastq`/`mkfastq` *before* Cell Ranger; cell-level via the R1 barcode, done
*by* Cell Ranger. The poly-dT is never sequenced -- R1 stops at exactly 16+12.

**SRA has no enforced convention.** `--split-files` writes `_1`,`_2`,`_3`. Two
files usually means R1,R2; three usually I1,R1,R2 -- but ordering varies by
submitter. **Verify by read length, never by numbering.**

**The classic trap:** SRA often flags the barcode read as *technical*, and
`fasterq-dump` skips technical reads by default. You get one file, assume the
deposit is broken. You need `--include-technical --split-files`.

Some submitters deposit Cell Ranger BAMs (recoverable with `bamtofastq`); a few
strip the barcode read entirely, and that data cannot be reprocessed at all.

**Cell Ranger's rigid naming** (`Sample_S1_L001_R1_001.fastq.gz`) is metadata, not
pedantry: it encodes sample name, sample index, **lane**, and read type. Lanes
matter because one library is often sequenced across several and must be merged;
`_001` handles size-split chunks. Point at a directory of many samples across
many lanes and it sorts itself out.

## 9. The two risky flags

```bash
cellranger count --id=SRR14575500 --transcriptome=refdata-gex-GRCh38-2024-A \
  --fastqs=fastqs/ --sample=SRR14575500 \
  --expect-cells=5000 --chemistry=SC3Pv3 --create-bam=false
```

**`--expect-cells` is nearly harmless now.** In v2 it drove the threshold
directly. Since v3, EmptyDrops estimates from the data and this is only a soft
prior -- pass 5,000 on a 15,000-cell run and you still get ~15,000. The dangerous
sibling is **`--force-cells`**, which bypasses cell calling entirely.

**`--chemistry` wrong is either loud or invisible.** It governs only how R1 is
sliced; alignment is untouched, since R2 is parsed independently. Declaring v3
(16+12) on a v2 library (16+10):

- *Usually aborts* -- R1 is 26 bp, 28 expected, and the whitelists differ so the
  valid-barcode fraction collapses. Good outcome.
- *If lengths permit it through* -- barcode (bases 1-16) is still fine, but the
  UMI becomes the true 10 bp plus 2 junk bases. The same molecule now yields
  different 12-mers, deduplication stops collapsing them, and **counts inflate.**

Modern Cell Ranger defaults to `--chemistry=auto` and detects by whitelist match
rate. Hardcoding it, as the tutorial does, overrides a safety net.

**Pipeline response:** `TENX_READ_DETECT` detects chemistry from R1 length, and
`QC_ASSERT` checks the observed valid-barcode fraction rather than trusting the
declaration.

## 10. Outputs

`outs/` holds `web_summary.html`, `metrics_summary.csv`,
`filtered_feature_bc_matrix{,.h5}`, `raw_feature_bc_matrix{,.h5}`.

The MTX triplet:

```
%%MatrixMarket matrix coordinate integer general
32285 9726 24857392      <- genes  cells  nonzeros
1247  3    5             <- gene 1247, cell 3, count 5
```

Sparsity is real -- a cell expresses 1,000-5,000 of ~30,000 genes, so **95-98% of
entries are zero** and a dense CSV would be gigabytes of nothing.

But **three files is convention, not engineering.** Matrix Market is a mid-90s
NIST text format with no metadata capability, so labels went into sidecars
(`barcodes.tsv.gz` in column order, `features.tsv.gz` in row order with
gene_id / gene_symbol / feature_type). Bioconductor's `Matrix` package used it,
10x shipped against that tooling in 2016, and it stuck. Cell Ranger already emits
the same data as a single `.h5` -- strictly better.

**So this pipeline standardises on HDF5 internally (`.h5ad`)** and reads MTX only
for compatibility.

Traps to design around:
- **Gene symbols are not unique.** Key on `gene_id`; call `var_names_make_unique()`.
- `feature_type` is how CITE-seq antibody counts and hashtags ride in the same matrix.

## 11. FastQC on R2 only

FastQC on R1 does not crash -- it *lies*. Its model assumes reads are independent
draws from a library, and R1 violates that by construction:

- **Sequence Duplication Levels** goes red: every read from one cell shares a
  barcode. That is the assay working.
- **Overrepresented sequences** flags the commonest barcodes. Also meaningless.
- **Adapter/k-mer content** is noise on a 28 bp synthetic construct.

Barcode quality lives in Cell Ranger's own metrics (`Q30 Bases in Barcode`,
`Fraction Valid Barcodes`). On R2, FastQC earns its place: adapter read-through,
end-of-read quality decay, poly-G from NovaSeq two-colour chemistry, rRNA
over-representation.

---

## Design decisions taken

| Decision | Choice | Why |
|---|---|---|
| Orchestration | Nextflow DSL2, thin processes over Python/R scripts | logic stays runnable standalone |
| Quantifiers | Cell Ranger + alevin-fry + kallisto\|bustools | one alignment-based, two lightweight; STARsolo dropped as redundant with Cell Ranger |
| Interface | all quantifiers emit `[meta, filtered, raw, metrics]` | adding one later is a module file, not a refactor |
| Chemistry | auto-detected, samplesheet overrides, mismatch is an error | the silent failure mode inflates counts |
| Part 1 boundary | closes with `COUNTS_TO_H5AD` | Part 2 never learns which quantifier ran |
| Cell Ranger | never containerised; bind-mounted from `params.cellranger_path` | 10x licence forbids redistribution |
| Container | one image, multi-arch, grown per lecture | avoids a 10 GB build for tools used at Part 12 |
| Samplesheet | experiment manifest, reserved columns + arbitrary passthrough | design declared once; no re-annotation at Part 3/5 |
