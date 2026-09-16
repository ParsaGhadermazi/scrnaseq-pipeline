#!/usr/bin/env python3
"""Fail a run on a bad quantification instead of passing a garbage matrix on.

The specific failure this catches: a wrong chemistry declaration mis-slices R1,
the barcodes stop matching the whitelist, and the valid-barcode fraction
collapses. Cell Ranger will happily continue and emit a near-empty or inflated
matrix. See docs/lectures/01-fastq-to-counts.md section 9.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

__version__ = "0.1.0"


def _num(x) -> float | None:
    try:
        s = str(x).strip().replace(",", "")
        return float(s[:-1]) / 100 if s.endswith("%") else float(s)
    except (TypeError, ValueError):
        return None


def parse_metrics(d: Path, quantifier: str) -> dict:
    """Pull a comparable set of numbers out of each tool's own report format.

    `d` is a metrics DIRECTORY -- every quantifier module stages its reports
    into one, so this contract does not change when a quantifier is added.
    Missing values come back as None and are reported "unavailable" rather than
    silently passing the check.
    """
    def _json(name):
        f = d / name
        return json.loads(f.read_text()) if f.exists() else {}

    if quantifier == "cellranger":
        f = d / "metrics_summary.csv"
        if not f.exists():
            return {}
        with open(f) as fh:
            row = next(csv.DictReader(fh))
        return {
            "valid_barcodes":        _num(row.get("Valid Barcodes")),
            "estimated_cells":       _num(row.get("Estimated Number of Cells")),
            "median_genes_per_cell": _num(row.get("Median Genes per Cell")),
            "sequencing_saturation": _num(row.get("Sequencing Saturation")),
            "reads_in_cells":        _num(row.get("Fraction Reads in Cells")),
            "mapped":                _num(row.get("Reads Mapped Confidently to Transcriptome")),
        }

    if quantifier == "kallisto":
        # inspect.json holds the barcode stats; run_info.json the alignment rate.
        ins, run = _json("inspect.json"), _json("run_info.json")
        pct = _num(ins.get("percentageReadsOnOnlist"))
        return {
            "valid_barcodes":  pct / 100 if pct is not None else None,
            "estimated_cells": _num(ins.get("numBarcodesOnOnlist")),
            "median_umis_per_cell": _num(ins.get("medianUMIsPerBarcode")),
            "mapped": (lambda x: x / 100 if x is not None else None)(_num(run.get("p_pseudoaligned"))),
        }

    # alevin-fry: simpleaf writes several small json reports
    log = _json("simpleaf_quant_log.json")
    gpl = _json("generate_permit_list.json")
    qnt = _json("quant.json")
    total = _num(gpl.get("total_reads"))
    corrected = _num(gpl.get("num_reads_with_corrected_barcode"))
    return {
        "valid_barcodes":  (corrected / total) if (total and corrected is not None) else None,
        "estimated_cells": _num(qnt.get("num_quantified_cells")) or _num(gpl.get("num_cells")),
        "mapped":          _num(log.get("percent_mapped")),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metrics", required=True, type=Path,
                    help="directory of the quantifier's own report files")
    ap.add_argument("--readinfo", required=True, type=Path)
    ap.add_argument("--quantifier", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--min-valid-barcodes", type=float, default=0.70)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--strict", action="store_true", help="exit non-zero on failure")
    ap.add_argument("--version", action="version", version=__version__)
    args = ap.parse_args()

    metrics = parse_metrics(args.metrics, args.quantifier)
    readinfo = json.loads(args.readinfo.read_text())

    checks, failures = {}, []

    vb = metrics.get("valid_barcodes")
    if vb is None:
        checks["valid_barcodes"] = {"status": "unavailable"}
    else:
        ok = vb >= args.min_valid_barcodes
        checks["valid_barcodes"] = {"value": vb, "threshold": args.min_valid_barcodes,
                                    "status": "pass" if ok else "FAIL"}
        if not ok:
            failures.append(
                f"valid barcode fraction {vb:.1%} < {args.min_valid_barcodes:.0%}. "
                f"Chemistry used was {readinfo.get('chemistry')} "
                f"(R1 = {readinfo.get('r1_length')} bp). A collapsed valid-barcode "
                "fraction is the signature of a wrong chemistry or a barcode read "
                "that is not what it claims to be."
            )

    if readinfo.get("chemistry_declared") not in (None, "auto") and \
       readinfo["chemistry_declared"] != readinfo["chemistry_detected"]:
        failures.append(f"declared {readinfo['chemistry_declared']} but R1 length implies "
                        f"{readinfo['chemistry_detected']}")

    report = {"sample": args.sample, "quantifier": args.quantifier,
              "metrics": metrics, "readinfo": readinfo,
              "checks": checks, "failures": failures,
              "status": "FAIL" if failures else "PASS",
              "qc_assert_version": __version__}
    args.out.write_text(json.dumps(report, indent=2))

    for f in failures:
        print(f"QC FAIL [{args.sample}]: {f}", file=sys.stderr)
    if failures and args.strict:
        sys.exit(1)
    print(f"QC {report['status']} [{args.sample}]")


if __name__ == "__main__":
    main()
