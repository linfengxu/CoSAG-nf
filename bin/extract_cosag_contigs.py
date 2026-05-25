#!/usr/bin/env python3
"""
Extract CoSAG contigs for GTDB input from cluster_data_tnf.json.

Selection rules:
1) If cluster has co_assembly_results.tnf_optimized.checkm2_optimized:
   use its completeness/contamination/contigs_file.
2) Otherwise use checkm2_coassembly + coassembly_contigs.
3) Keep only contigs passing:
   completeness >= min_completeness AND contamination <= max_contamination.

Outputs:
- One renamed FASTA per passing cluster: cluster_<cluster_id>.fasta
- manifest TSV with source and QC fields
"""

import argparse
import json
import shutil
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


def _to_float(v: Any) -> Optional[float]:
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def choose_contig(cluster: Dict[str, Any]) -> Tuple[Optional[str], Optional[float], Optional[float], Optional[str]]:
    """
    Return (contig_path, completeness, contamination, source).
    source in {"optimized", "coassembly"}.
    """
    coasm = cluster.get("co_assembly_results", {}) or {}
    tnf = coasm.get("tnf_optimized", {}) or {}
    opt = tnf.get("checkm2_optimized", {}) or {}

    opt_contig = opt.get("contigs_file")
    opt_comp = _to_float(opt.get("completeness"))
    opt_contam = _to_float(opt.get("contamination"))

    if opt_contig and opt_comp is not None and opt_contam is not None:
        return str(opt_contig), opt_comp, opt_contam, "optimized"

    base = cluster.get("checkm2_coassembly", {}) or {}
    base_contig = cluster.get("coassembly_contigs")
    base_comp = _to_float(base.get("completeness"))
    base_contam = _to_float(base.get("contamination"))

    if base_contig and base_comp is not None and base_contam is not None:
        return str(base_contig), base_comp, base_contam, "coassembly"

    return None, None, None, None


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--json", required=True, type=Path, help="cluster_data_tnf.json path")
    p.add_argument("--outdir", required=True, type=Path, help="output folder for renamed FASTA files")
    p.add_argument("--min_completeness", type=float, default=50.0)
    p.add_argument("--max_contamination", type=float, default=10.0)
    p.add_argument(
        "--mode",
        choices=["copy", "symlink"],
        default="copy",
        help="How to place selected FASTA files into outdir",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()

    if not args.json.exists():
        print(f"ERROR: json not found: {args.json}", file=sys.stderr)
        return 1

    args.outdir.mkdir(parents=True, exist_ok=True)

    data = json.loads(args.json.read_text(encoding="utf-8"))
    clusters: List[Dict[str, Any]] = data.get("clusters", [])

    manifest_lines = [
        "cluster_id\tsource\tcompleteness\tcontamination\tinput_contig\toutput_contig\tselected"
    ]

    selected_n = 0
    skipped_n = 0

    for c in clusters:
        cid = str(c.get("cluster_id", "UNKNOWN"))
        contig, comp, contam, source = choose_contig(c)

        if not contig or comp is None or contam is None:
            manifest_lines.append(
                f"{cid}\tNA\tNA\tNA\tNA\tNA\tno_metrics_or_contig"
            )
            skipped_n += 1
            continue

        src = Path(contig)
        if not src.exists():
            manifest_lines.append(
                f"{cid}\t{source}\t{comp:.2f}\t{contam:.2f}\t{contig}\tNA\tmissing_input_fasta"
            )
            skipped_n += 1
            continue

        passes = comp >= args.min_completeness and contam <= args.max_contamination
        out_name = f"cluster_{cid}.fasta"
        dst = args.outdir / out_name

        if passes:
            if args.mode == "copy":
                shutil.copy2(src, dst)
            else:
                if dst.exists() or dst.is_symlink():
                    dst.unlink()
                dst.symlink_to(src.resolve())
            manifest_lines.append(
                f"{cid}\t{source}\t{comp:.2f}\t{contam:.2f}\t{contig}\t{dst}\tpass"
            )
            selected_n += 1
        else:
            manifest_lines.append(
                f"{cid}\t{source}\t{comp:.2f}\t{contam:.2f}\t{contig}\tNA\tfail_threshold"
            )
            skipped_n += 1

    manifest = args.outdir / "selected_contigs.tsv"
    manifest.write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")

    print(f"Selected: {selected_n}")
    print(f"Skipped : {skipped_n}")
    print(f"Manifest: {manifest}")

    if selected_n == 0:
        print("ERROR: no contigs selected", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
