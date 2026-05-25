#!/usr/bin/env python3
"""
update_cluster_optimized.py
----------------------------
Update cluster_data_updated.json with TNF optimization results
from per-cluster *_tnf_optimized.json files.

Adds co_assembly_results.tnf_optimized section to each cluster.
"""

import argparse
import json
import logging
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def load_json(path: Path) -> Dict:
    with path.open(encoding="utf-8") as fh:
        return json.load(fh)


def find_cluster(data: Dict, cluster_id: str) -> Optional[Dict]:
    for c in data.get("clusters", []):
        if str(c.get("cluster_id")) == str(cluster_id):
            return c
    return None


def extract_cluster_id(path: Path, opt: Optional[Dict] = None) -> Optional[str]:
    """
    Extract cluster ID with priority:
      1) opt['cluster_id'] (most reliable)
      2) filename like cluster_R1_C408_optimized.json
    """
    if isinstance(opt, dict) and opt.get("cluster_id"):
        return str(opt["cluster_id"])

    stem = path.stem
    m = re.match(r"^cluster_(.+?)_optimized(?:_\d+)?$", stem)
    if m:
        return m.group(1)
    return None


# ---------------------------------------------------------------------------
# Core update
# ---------------------------------------------------------------------------
def update_with_tnf(cluster: Dict, opt: Dict) -> None:
    cluster.setdefault("co_assembly_results", {})
    coasm = cluster["co_assembly_results"]

    # Keep original co-assembly CheckM2 metrics untouched for comparison.
    if "original_checkm2_coassembly" not in coasm and isinstance(cluster.get("checkm2_coassembly"), dict):
        coasm["original_checkm2_coassembly"] = dict(cluster["checkm2_coassembly"])

    # Compatible with both "full TNF summary" JSON and current optimizer output JSON.
    selected_sags = opt.get("selected_sags")
    if not isinstance(selected_sags, list) or not selected_sags:
        members = opt.get("members", opt.get("sags", []))
        selected_sags = members if isinstance(members, list) else []

    selected_sag_count = opt.get("selected_sag_count")
    if selected_sag_count is None:
        selected_sag_count = opt.get("member_count", len(selected_sags))

    initial_sag_count = opt.get("initial_sag_count")
    if initial_sag_count is None:
        leaf_sags = cluster.get("leaf_sags", [])
        if isinstance(leaf_sags, list) and leaf_sags:
            initial_sag_count = len(leaf_sags)
        else:
            initial_sag_count = cluster.get("cluster_size", selected_sag_count)

    optimization_summary = opt.get("optimization_summary", {})
    converged = opt.get("converged")
    if converged is None and isinstance(optimization_summary, dict):
        converged = optimization_summary.get("meets_targets")

    tnf_section = {
        "cluster_id":            opt.get("cluster_id", cluster.get("cluster_id")),
        "selected_sags":         selected_sags,
        "selected_sag_count":    selected_sag_count,
        "initial_sag_count":     initial_sag_count,
        "outliers_removed":      opt.get("outliers_removed", []),
        "iterative_removed":     opt.get("iterative_removed", []),
        "final_mean_tnf_dist":   opt.get("final_mean_tnf_dist"),
        "converged":             converged,
        "distance_history":      opt.get("distance_history", []),
        "removal_log":           opt.get("removal_log", []),
        "cfg":                   opt.get("cfg", {}),
        "optimization_summary":  optimization_summary if isinstance(optimization_summary, dict) else {},
        "checkm2_optimized": {
            "completeness":      opt.get("completeness"),
            "contamination":     opt.get("contamination"),
            "genome_size":       opt.get("genome_size"),
            "contig_n50":        opt.get("contig_n50"),
            "gc_content":        opt.get("gc_content"),
            "coding_density":    opt.get("coding_density"),
            "contigs_file":      opt.get("contigs_file"),
            "result_type":       opt.get("result_type"),
            "iteration":         opt.get("iteration"),
            "member_count":      opt.get("member_count"),
            "optimization_tag":  opt.get("tag"),
        },
    }
    coasm["tnf_optimized"] = tnf_section


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def run(
    main_json:  Path,
    opt_jsons:  List[Path],
    out_json:   Path,
) -> None:

    data = load_json(main_json)
    updated = []

    for opt_path in opt_jsons:
        if not opt_path.exists():
            log.warning("File not found: %s", opt_path)
            continue

        opt = load_json(opt_path)
        cid = extract_cluster_id(opt_path, opt)
        if cid is None:
            log.warning("Cannot extract cluster ID from: %s", opt_path.name)
            continue

        cluster = find_cluster(data, cid)
        if cluster is None:
            log.warning("Cluster %s not found in main JSON", cid)
            continue

        update_with_tnf(cluster, opt)
        updated.append(cid)
        selected_count = opt.get("selected_sag_count")
        if selected_count is None:
            selected_count = opt.get("member_count", len(opt.get("selected_sags", opt.get("members", [])) or []))
        initial_count = opt.get("initial_sag_count")
        if initial_count is None:
            initial_count = len(cluster.get("leaf_sags", [])) if isinstance(cluster.get("leaf_sags", []), list) else cluster.get("cluster_size", 0)
        converged = opt.get("converged")
        if converged is None and isinstance(opt.get("optimization_summary"), dict):
            converged = opt["optimization_summary"].get("meets_targets")
        log.info(
            "Cluster %s: selected=%d/%d  converged=%s  final_dist=%.4f",
            cid,
            selected_count,
            initial_count,
            converged,
            opt.get("final_mean_tnf_dist", 0),
        )

    out_json.write_text(json.dumps(data, indent=2, ensure_ascii=False))
    log.info("Updated %d clusters → %s", len(updated), out_json)

    if not updated:
        log.warning("No clusters updated; writing unchanged JSON")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--main_json",   required=True, type=Path,
                   help="cluster_data_updated.json")
    p.add_argument("--opt_jsons",   required=True, nargs="+", type=Path,
                   help="cluster_*_tnf_optimized.json files")
    p.add_argument("--out_json",    required=True, type=Path)
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(args.main_json, args.opt_jsons, args.out_json)
