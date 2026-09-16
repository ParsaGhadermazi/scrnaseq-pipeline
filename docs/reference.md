# References and index building

**Rule: if you do not supply an index, the pipeline builds one.**

| Quantifier | Supply this | ...or build from | Module |
|---|---|---|---|
| cellranger | `--cellranger_ref` (a `refdata-gex-*` dir) | `--genome_fasta` + `--gtf` | `CELLRANGER_MKREF` |
| alevinfry  | `--index` (simpleaf's output dir)          | `--genome_fasta` + `--gtf` | `SIMPLEAF_INDEX` |
| kallisto   | `--index` (dir with `index.idx`, `t2g.txt`) | `--genome_fasta` + `--gtf` | `KB_REF` |

## Caching

Indexes are **not** written to `--outdir`. They use Nextflow's `storeDir`, keyed
on `--index_cache` (default `<launchDir>/references`). If the index already
exists there, the build process is skipped entirely.

This matters because a human reference is expensive: `cellranger mkref` is hours
of CPU, needs ~32 GB+ RAM, and produces ~11 GB. Rebuilding per run would make
the pipeline unusable.

> **Point `--index_cache` at an external drive if your internal disk is tight.**

## Two choices that are easy to get wrong

### 1. GTF biotype filtering (`--gtf_biotypes`)

`cellranger mkgtf` runs before `mkref`, keeping only the biotypes listed in
`--gtf_biotypes` (10x's set by default) and stripping pseudogenes, readthrough
transcripts, and retained-intron variants.

This is **not cosmetic**. Cell Ranger counts a read only if it maps uniquely to
*one* gene; a read overlapping two annotated features is discarded as ambiguous.
So a pseudogene resembling its parent gene does not merely add a junk row -- it
**silently deletes counts from the real gene**. Feeding a raw Ensembl GTF is a
quiet accuracy bug, not a harmless shortcut.

### 2. Splici vs spliced-only (`--af_ref_type`, `--kb_workflow`)

Default is `spliced+intronic` (splici) for alevin-fry.

Separating spliced from unspliced reads is the **entire basis of RNA velocity**.
Cell Ranger merges them into one number and the distinction is gone forever --
which is why Part 13 (scVelo) cannot run off a plain Cell Ranger matrix.

Choosing a spliced-only index here permanently forecloses Part 13 for that data.
The equivalent for kallisto is `--kb_workflow nac`, which is **not** the default
(`standard` is), so set it if you want velocity on the kallisto path.

`--read_length` (default 91) must match your R2 length: it sets how much
intronic flank surrounds each exon-intron boundary in a splici index.

## Provenance

Every index build writes a `ref.json` recording the source FASTA, GTF, biotypes
kept, reference type, tool versions, and whether the index is velocity-capable.
Two runs against different references are not comparable, so this travels with
the counts rather than living only in a build log.

## Verification status

**None of the three index modules has been executed.** They are written against
each tool's documented interface but not run:

- `CELLRANGER_MKREF` -- needs the licensed x86_64 binary, ~32 GB+ RAM. Not
  runnable on Apple silicon or on a 16 GB laptop.
- `SIMPLEAF_INDEX` -- the splici path (`--fasta` + `--gtf`) is unrun. Only
  `--ref-seq` against a toy transcriptome has been exercised.
- `KB_REF` -- unrun. `KB_COUNT` was tested against a hand-built kallisto index.

Expect to debug these on first real use. See `docs/COVERAGE.md`.
