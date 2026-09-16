#!/usr/bin/env python3
"""Identify 10x read roles by observed length and rename to Cell Ranger convention.

Why this exists (docs/lectures/01-fastq-to-counts.md, sections 8-9):

  * SRA enforces no convention. `--split-files` writes _1/_2/_3 in an order that
    varies by submitter, so numbering cannot be trusted -- only read length can.
  * SRA frequently flags the barcode read as "technical"; without
    `--include-technical` fasterq-dump silently returns R2 only.
  * Declaring the wrong chemistry either aborts loudly or, if read lengths let
    it through, mis-slices the UMI so deduplication stops collapsing PCR
    duplicates and counts inflate invisibly.

Length -> role, for 3' gene expression:

    8, 10 bp   -> I1 / I2   sample index (already consumed by bcl2fastq)
    26 bp      -> R1        16 bp barcode + 10 bp UMI  => SC3Pv2
    28 bp      -> R1        16 bp barcode + 12 bp UMI  => SC3Pv3
    >= 50 bp   -> R2        cDNA
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import sys
from collections import Counter
from pathlib import Path

__version__ = "0.1.0"

N_PROBE = 2000          # records sampled per file to establish modal length
CHEMISTRY_BY_R1 = {26: "SC3Pv2", 28: "SC3Pv3"}
UMI_LEN = {"SC3Pv2": 10, "SC3Pv3": 12}
BARCODE_LEN = 16


def _open(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def modal_read_length(path: Path) -> tuple[int, int]:
    """Return (modal length, records sampled). FASTQ is 4 lines per record."""
    lengths: Counter[int] = Counter()
    with _open(path) as fh:
        for i, line in enumerate(fh):
            if i >= N_PROBE * 4:
                break
            if i % 4 == 1:
                lengths[len(line.rstrip("\n"))] += 1
    if not lengths:
        sys.exit(f"ERROR: {path} contains no reads")
    return lengths.most_common(1)[0][0], sum(lengths.values())


def classify(lengths: dict[str, int]) -> dict[str, str]:
    """Map filename -> role. Raises on anything ambiguous rather than guessing."""
    roles: dict[str, str] = {}
    for name, ln in lengths.items():
        if ln <= 12:
            roles[name] = "I1"
        elif ln in CHEMISTRY_BY_R1:
            roles[name] = "R1"
        elif ln >= 50:
            roles[name] = "R2"
        else:
            roles[name] = "UNKNOWN"

    # Second index read: disambiguate the two short ones by length then name.
    index_reads = sorted(n for n, r in roles.items() if r == "I1")
    for extra in index_reads[1:]:
        roles[extra] = "I2"

    n_r1 = sum(1 for r in roles.values() if r == "R1")
    n_r2 = sum(1 for r in roles.values() if r == "R2")

    if n_r1 == 0:
        raise SystemExit(
            "ERROR: no read of 26 or 28 bp found -- the cell barcode read is missing.\n"
            "  Most likely cause: fasterq-dump ran without --include-technical, so SRA\n"
            "  dropped the barcode read as 'technical'. Re-download with:\n"
            "      fasterq-dump --split-files --include-technical <SRR>\n"
            f"  Observed lengths: {lengths}"
        )
    if n_r1 > 1 or n_r2 > 1:
        raise SystemExit(
            f"ERROR: ambiguous read roles, cannot assign by length: {lengths}\n"
            "  Two files share a role. Supply fastq_1/fastq_2 explicitly in the samplesheet."
        )
    if n_r2 == 0:
        raise SystemExit(f"ERROR: no cDNA read (>=50 bp) found. Observed lengths: {lengths}")
    return roles


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fastqs", nargs="+", required=True, type=Path)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--json", required=True, type=Path)
    ap.add_argument("--declared-chemistry", default=None,
                    help="Samplesheet override. A mismatch is an error, not a warning.")
    ap.add_argument("--symlink", action="store_true",
                    help="Symlink instead of copy (saves disk; needs a stageable filesystem).")
    ap.add_argument("--version", action="version", version=__version__)
    args = ap.parse_args()

    lengths = {p.name: modal_read_length(p)[0] for p in args.fastqs}
    roles = classify(lengths)

    by_role = {r: n for n, r in roles.items()}
    r1_len = lengths[by_role["R1"]]
    detected = CHEMISTRY_BY_R1[r1_len]

    if args.declared_chemistry and args.declared_chemistry != "auto":
        if args.declared_chemistry != detected:
            raise SystemExit(
                f"ERROR: chemistry mismatch for sample '{args.sample}'.\n"
                f"  samplesheet declares : {args.declared_chemistry}\n"
                f"  R1 length implies    : {detected} (R1 = {r1_len} bp)\n"
                "  Declaring the wrong chemistry mis-slices the UMI: deduplication stops\n"
                "  collapsing PCR duplicates and counts inflate silently. Refusing to run.\n"
                "  Remove the 'chemistry' column to trust detection, or fix the declaration."
            )
        chemistry = args.declared_chemistry
    else:
        chemistry = detected

    # rename to the convention Cell Ranger parses as metadata:
    #   <sample>_S1_L001_<read>_001.fastq.gz
    args.outdir.mkdir(parents=True, exist_ok=True)
    renamed: dict[str, str] = {}
    for name, role in roles.items():
        if role == "UNKNOWN":
            continue
        src = next(p for p in args.fastqs if p.name == name)
        dst = args.outdir / f"{args.sample}_S1_L001_{role}_001.fastq.gz"
        if dst.exists():
            dst.unlink()
        if args.symlink:
            os.symlink(os.path.abspath(src), dst)
        elif str(src).endswith(".gz"):
            shutil.copy2(src, dst)
        else:
            with open(src, "rb") as fi, gzip.open(dst, "wb") as fo:
                shutil.copyfileobj(fi, fo)
        renamed[role] = dst.name

    info = {
        "sample": args.sample,
        "chemistry": chemistry,
        "chemistry_detected": detected,
        "chemistry_declared": args.declared_chemistry,
        "r1_length": r1_len,
        "barcode_length": BARCODE_LEN,
        "umi_length": UMI_LEN[detected],
        "observed_lengths": lengths,
        "roles": roles,
        "renamed": renamed,
        "detector_version": __version__,
    }
    args.json.write_text(json.dumps(info, indent=2))
    print(json.dumps(info, indent=2))


if __name__ == "__main__":
    main()
