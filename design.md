# Design principles

This document describes the design principles for this RNA-seq pipeline. It was
started as a short brief and has been filled in as decisions were actually made
and tested; every rule below was a deliberate choice, and the reason is recorded
because the reasons are what generalise.

## Top level

The user selects the data type:

- **Single-cell RNA-seq** — built (stages `counts`, `qc`)
- **Single-cell ATAC-seq** — planned, after the scRNA-seq track completes
- **Bulk RNA-seq** — planned

## Structural rules

**Nextflow DSL2 orchestrates; logic lives in scripts.** Processes are thin
wrappers around `bin/*.py` and `bin/*.R`, so every step is runnable standalone
and testable outside the workflow engine.

**Stages run one at a time.** `counts` → `qc` → `annotate` → `de`. Each stage is
a self-contained invocation that locates the previous stage's output by
convention. This exists because analysts iterate on one step at a time; putting
clustering inside `qc` would mean re-running ambient correction to change one
resolution parameter.

**One samplesheet is the experiment manifest**, declared once and reused by every
stage. Reserved columns are interpreted; anything else is carried verbatim into
the meta map and lands in `adata.obs`. Re-annotating samples at each stage is
where samples get mismatched.

**`.h5ad` is the interchange format**, always. Not MTX — three files is a 1990s
convention with no metadata capability, not an engineering choice.

**Interfaces before implementations.** Every quantifier emits its native output
and a single converter knows the per-tool layouts. Adding a tool is a module plus
a reader function, never a refactor.

## Tool selection

**Implement every tool the tutorial uses. Add a cross-language alternative only
when it is genuinely competitive — never for symmetry.** If the best tool for a
step exists only in R, use R and say so. This is why Part 1 is Python-led and
Part 2 is R-led. Tracked per step in `docs/COVERAGE.md`.

**One aligner per run; one ambient method per run; doublet callers may be a
list.** The rule is cost: ambient correction rewrites the entire matrix, so a
second method doubles storage. Doublet callers emit a per-cell CSV, so running
several is nearly free — which is what makes consensus removal practical.

## Correctness rules

**Fail loudly rather than emit quietly-wrong output.** A collapsed valid-barcode
fraction, a chemistry declaration contradicting the read length, a `batch_key`
perfectly confounded with `condition_key` — all stop the run. The failure mode
that matters is not the crash, it is the plausible-looking wrong answer.

**Record what was actually used, not what was requested.** Versions are captured
from running binaries, the container digest is recorded, and the environment lock
file is generated *from the built image*. A pin that silently does not bind is
worse than no pin.

**Never guess availability from a failed probe.** A network timeout and "no build
for this architecture" are different things; conflating them silently changed
what went into the image.

**Raw counts are never overwritten.** `layers["counts"]` survives every stage,
because differential expression must never run on ambient-corrected or
batch-integrated values.

## Constraints that shaped the build

**Cell Ranger is never containerised** — 10x's licence forbids redistribution. It
is bind-mounted from a user-supplied path.

**R can read `.h5ad` but cannot write it.** No R module writes `.h5ad`; the one
step that returns a matrix writes MTX for a Python step to convert. See
`docs/python-r-bridge.md` — this was established by testing, after two
plausible-looking approaches turned out not to exist.
