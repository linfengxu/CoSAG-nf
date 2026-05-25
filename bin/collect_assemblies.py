#!/usr/bin/env python3
"""
collect_assemblies.py
---------------------
Collect optimized CoSAG contigs and co-assembly contigs,
apply quality filtering, and generate a final collection report.

Input sources:
    --optimized_jsons : cluster_*_optimized.json from COSAG_OPTIMIZATION
    --coassembly_json : cluster_data_updated.json from UPDATE_CLUSTER_JSON
    --out_dir         : output directory

Output directories:
    bset_cosag/          : TNF-optimized best contigs (OPTIMAL / FALLBACK tier)
    co_assembly/         : direct co-assembly contigs passing QC
    all_accepted_cosags/ : union of above (for dRep / GTDB)
    filtered_out/        : contigs failing QC
"""

import argparse
import json
import shutil
import sys
from pathlib import Path


DEFAULTS = {
    "min_completeness": 50.0,   # for co_assembly contigs (not TNF-optimized)
    "max_contamination": 10.0,
}


def load_json(path: Path) -> dict:
    with path.open() as fh:
        return json.load(fh)


def collect_optimized(
    json_files: list,
    bset_dir: Path,
    all_dir: Path,
    filtered_dir: Path,
    cfg: dict,
) -> tuple[list, list]:
    """Collect TNF-optimized contigs. OPTIMAL and FALLBACK tiers are accepted."""
    passed, failed = [], []

    for jf in json_files:
        p = Path(jf.strip())
        if not p.exists():
            continue
        data = load_json(p)
        tier    = data.get("tier", "BEST_AVAILABLE")
        cid     = data.get("cluster_id", "unknown")
        comp    = data.get("completeness", 0)
        cont    = data.get("contamination", 100)
        gsize   = data.get("genome_size", 0)
        n50     = data.get("contig_n50", 0)
        contigs = Path(data.get("contigs", ""))

        if not contigs.exists():
            # try all_contig directory
            parent = p.parent / "all_contig" / \
                     f"cluster_{cid}_best_contigs.fasta"
            if parent.exists():
                contigs = parent
            else:
                print(f"WARNING: contigs not found for cluster {cid}", file=sys.stderr)
                continue

        accepted = tier in ("OPTIMAL", "FALLBACK")
        dst_name = f"cluster_{cid}_CoSAG_optimized.fasta"

        info = {
            "cluster_id": cid,
            "tier":        tier,
            "completeness": comp,
            "contamination": cont,
            "genome_size":  gsize,
            "contig_n50":   n50,
            "source":       "tnf_optimized",
            "filename":     dst_name,
        }

        if accepted:
            shutil.copy2(contigs, bset_dir / dst_name)
            shutil.copy2(contigs, all_dir  / dst_name)
            passed.append(info)
        else:
            shutil.copy2(contigs, filtered_dir / dst_name)
            info["reason"] = f"tier={tier}"
            failed.append(info)

    return passed, failed


def collect_coassembly(
    updated_json: Path,
    coasm_dir: Path,
    all_dir: Path,
    filtered_dir: Path,
    cfg: dict,
) -> tuple[list, list]:
    """Collect direct co-assembly contigs that pass QC thresholds."""
    if not updated_json or not updated_json.exists():
        return [], []

    data    = load_json(updated_json)
    passed, failed = [], []

    for cluster in data.get("clusters", []):
        cid     = cluster["cluster_id"]
        r1      = cluster.get("co_assembly_results", {}).get("round_1", {})
        comp    = r1.get("completeness",  0)
        cont    = r1.get("contamination", 100)
        gsize   = r1.get("genome_size",   0)
        n50     = r1.get("contig_n50",    0)
        src     = r1.get("assembly_file", "")

        if not src or not Path(src).exists():
            continue

        dst_name = f"cluster_{cid}_CoSAG.fasta"
        reasons  = []

        if cont > cfg["max_contamination"]:
            reasons.append(f"contamination {cont:.2f}% > {cfg['max_contamination']}%")
        if comp < cfg["min_completeness"]:
            reasons.append(f"completeness {comp:.2f}% < {cfg['min_completeness']}%")

        info = {
            "cluster_id":   cid,
            "completeness": comp,
            "contamination": cont,
            "genome_size":  gsize,
            "contig_n50":   n50,
            "source":       "co_assembly",
            "filename":     dst_name,
        }

        if not reasons:
            shutil.copy2(src, coasm_dir / dst_name)
            shutil.copy2(src, all_dir   / dst_name)
            passed.append(info)
        else:
            shutil.copy2(src, filtered_dir / dst_name)
            info["reason"] = "; ".join(reasons)
            failed.append(info)

    return passed, failed


def write_report(
    out_path: Path,
    optimized_pass: list, optimized_fail: list,
    coasm_pass: list,     coasm_fail: list,
    cfg: dict,
) -> None:
    total_pass = len(optimized_pass) + len(coasm_pass)
    total_fail = len(optimized_fail) + len(coasm_fail)

    with out_path.open("w") as fh:
        fh.write("=" * 80 + "\n")
        fh.write("CoSAG Assembly Collection Report\n")
        fh.write("=" * 80 + "\n\n")
        fh.write("QC thresholds (co-assembly):\n")
        fh.write(f"  min_completeness:  {cfg['min_completeness']}%\n")
        fh.write(f"  max_contamination: {cfg['max_contamination']}%\n\n")
        fh.write(f"TNF-optimized (OPTIMAL/FALLBACK):  {len(optimized_pass)} passed, "
                 f"{len(optimized_fail)} failed\n")
        fh.write(f"Co-assembly direct:                {len(coasm_pass)} passed, "
                 f"{len(coasm_fail)} failed\n")
        fh.write(f"Total in all_accepted_cosags/:     {total_pass}\n")
        fh.write(f"Total filtered_out/:               {total_fail}\n\n")

        for section, items in [
            ("TNF-Optimized Passed", optimized_pass),
            ("Co-Assembly Passed",   coasm_pass),
            ("TNF-Optimized Failed", optimized_fail),
            ("Co-Assembly Failed",   coasm_fail),
        ]:
            if not items:
                continue
            fh.write("-" * 60 + "\n")
            fh.write(f"{section} ({len(items)})\n")
            fh.write("-" * 60 + "\n")
            for i in sorted(items, key=lambda x: -x["completeness"]):
                fh.write(f"  {i['filename']}\n")
                fh.write(f"    comp={i['completeness']:.1f}%  "
                         f"cont={i['contamination']:.1f}%  "
                         f"size={i['genome_size']:,}bp  "
                         f"n50={i['contig_n50']:,}bp\n")
                if "reason" in i:
                    fh.write(f"    reason: {i['reason']}\n")
                if "tier" in i:
                    fh.write(f"    tier: {i['tier']}\n")
            fh.write("\n")

    print(f"Report: {out_path}")
    print(f"Total accepted: {total_pass}  |  Filtered: {total_fail}")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--optimized_jsons", nargs="*", default=[],
                   help="cluster_*_optimized.json files from COSAG_OPTIMIZATION")
    p.add_argument("--coassembly_json", type=Path, default=None,
                   help="cluster_data_updated.json from UPDATE_CLUSTER_JSON")
    p.add_argument("--out_dir",         required=True, type=Path)
    p.add_argument("--min_completeness", type=float,
                   default=DEFAULTS["min_completeness"])
    p.add_argument("--max_contamination", type=float,
                   default=DEFAULTS["max_contamination"])
    return p.parse_args()


def main():
    args = parse_args()
    cfg  = {k: getattr(args, k) for k in DEFAULTS}

    bset_dir     = args.out_dir / "bset_cosag"
    coasm_dir    = args.out_dir / "co_assembly"
    all_dir      = args.out_dir / "all_accepted_cosags"
    filtered_dir = args.out_dir / "filtered_out"

    for d in (bset_dir, coasm_dir, all_dir, filtered_dir):
        d.mkdir(parents=True, exist_ok=True)

    opt_pass, opt_fail = collect_optimized(
        args.optimized_jsons, bset_dir, all_dir, filtered_dir, cfg)

    ca_pass, ca_fail = collect_coassembly(
        args.coassembly_json, coasm_dir, all_dir, filtered_dir, cfg)

    write_report(
        args.out_dir / "collection_summary.txt",
        opt_pass, opt_fail, ca_pass, ca_fail, cfg,
    )

    if len(opt_pass) + len(ca_pass) == 0:
        print("ERROR: No contigs passed QC!", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()