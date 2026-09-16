# The Python/R bridge

The pipeline's interchange format is `.h5ad`. Python-side steps read and write it
natively via `anndata`. The R side is **asymmetric**, and this was established by
testing rather than by reading documentation -- two plausible-looking approaches
turned out to be fiction.

## What actually works

| Direction | Mechanism | Works offline? |
|---|---|---|
| R **reads** `.h5ad` | `zellkonverter::readH5AD(file, reader = "R")` | **yes** -- native rhdf5 |
| R **writes** `.h5ad` | *no native path exists* | **no** |

`reader` is a real argument of `readH5AD` (default `c("python", "R")`). Passing
`reader = "R"` uses rhdf5 and never touches Python. Verified with
`docker run --network none`.

## The trap

`writeH5AD()` has **no** `writer` argument. Its formals are:

    sce, file, X_name, skip_assays, compression, version, ...

Passing `writer = "R"` is silently absorbed by `...`, and the function runs its
**default Python path**. That bootstraps basilisk/reticulate, which clones pyenv
and **compiles CPython from source inside the container at runtime**. It is
non-reproducible, enormous, and fails outright with no network:

    writeH5AD ... -> pyenv_bootstrap -> download.file
    cannot open URL '.../pyenv-installer' : Could not resolve hostname

There is also no environment variable that disables this. Setting
`ZELLKONVERTER_USE_BASILISK=0` does nothing -- that variable does not exist.

## The rule

> **R modules READ `.h5ad` with `reader = "R"`. No R module WRITES `.h5ad`.**

Any R step that must return a **matrix** writes Matrix Market plus TSVs, and a
small Python step converts it to `.h5ad`:

```r
writeMM(counts, "out/matrix.mtx")
write.table(colnames(sce), "out/barcodes.tsv", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(sce), "out/features.tsv", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.csv(as.data.frame(colData(sce)), "out/obs.csv")
```

This costs one conversion in the entire DAG, because **only ambient correction
rewrites a matrix**. Every other R step (`emptyDrops`, `scDblFinder`,
`DoubletFinder`) emits per-cell annotations as CSV, which Python merges into
`obs` -- no matrix ever crosses the boundary.

Verified round-trip: R read a 50x20 `.h5ad`, doubled the counts, wrote MTX;
Python rebuilt it with matching shape, exactly 2x the original sum, and `obs`
columns preserved -- all with networking disabled.

## Alternative not pursued

`anndataR` (Bioconductor) offers native R read *and* write of `.h5ad`, and
zellkonverter's own startup banner now recommends it. It was not adopted because
the MTX route is already verified, needs no new dependency, and is exercised
exactly once per run. Worth revisiting if R-side matrix writing ever becomes
common.


---

# Upstream version traps

Not every failure is in our code. These are incompatibilities between installed
library versions, found by measurement and pinned here so they are not
rediscovered.

## scanpy's `harmony_integrate` is stale relative to harmonypy 2.0.0

`sc.external.pp.harmony_integrate` ends with:

```python
adata.obsm[adjusted_basis] = harmony_out.Z_corr.T
```

That transpose was correct for **harmonypy 1.x**, where `Z_corr` was
`(n_pcs, n_cells)`. **harmonypy 2.0.0 returns `Z_corr` as `(n_cells, n_pcs)`**,
so the wrapper transposes a correctly-oriented matrix into the wrong shape and
AnnData rejects it:

```
ValueError: Value passed for key 'X_integrated' is of incorrect shape.
Values of obsm must match dimensions ('obs',) of parent.
Value had shape (50,) while it should have had (1930,).
```

Measured directly: input `X_pca` `(1930, 50)` -> `Z_corr` `(1930, 50)` ->
`Z_corr.T` `(50, 1930)`.

**What we do:** `bin/integrate.py` calls `harmonypy.run_harmony` directly and
picks the orientation by checking which axis equals `n_obs`, rather than
trusting either version's convention:

```python
Z = np.asarray(ho.Z_corr)
emb = Z if Z.shape[0] == adata.n_obs else Z.T if Z.shape[1] == adata.n_obs else error
```

This is version-agnostic: it works whichever way a future harmonypy orients its
output, and fails loudly if neither axis matches.

Verified on real data: `(1930, 50)` embedding, distinct from `X_pca`, batch
silhouette -0.0079 -> -0.0139.


## Cross-architecture resolution

The image is published multi-arch (`linux/amd64` + `linux/arm64`) from the same
`environment.*.yml` pins. Both halves resolved to **identical package versions**
-- verified by reading `environment.lock.yml` out of each published image:

| package | linux/amd64 | linux/arm64 |
|---|---|---|
| sra-tools | 3.4.1 | 3.4.1 |
| salmon | 2.7.0 | 2.7.0 |
| alevin-fry | 0.18.3 | 0.18.3 |
| simpleaf | 0.30.0 | 0.30.0 |
| kallisto | 0.52.0 | 0.52.0 |
| bustools | 0.45.1 | 0.45.1 |
| samtools | 1.24 | 1.24 |
| r-base | 4.5.3 | 4.5.3 |
| scanpy / anndata | 1.11.5 / 0.12.19 | 1.11.5 / 0.12.19 |

Only the conda **build strings** differ (`hd612981_0` vs `h6789a04_0`), which is
correct -- they encode the toolchain and architecture, not the software version.

`environment.lock.yml` is generated *inside* each image for exactly this reason:
the lock is necessarily per-architecture, and it is the artifact that actually
pins a run. The pipeline version, container digest and per-tool versions in the
run manifest sit on top of it.

**Versions do drift across time, though, not across architecture.** An arm64
image built weeks earlier carried salmon 2.6.0, alevin-fry 0.18.0 and simpleaf
0.28.0; this build picked up 2.7.0, 0.18.3 and 0.30.0 from the same `>=` pins.
That is the argument for pinning against a published lock file rather than
re-solving, if byte-identical reruns ever matter.

**On the architecture worry that shaped early decisions:** bioconda's
`linux-aarch64` coverage turned out to be as complete as `linux-64` for this
stack. The prediction that `sra-tools`, `kallisto` and `salmon` would need source
builds on arm64 was wrong -- all resolve from conda on both. The only two source
builds are properties of the packages, not the architecture: `kb-python` is a pip
wrapper, and `DoubletFinder` is GitHub-only.
