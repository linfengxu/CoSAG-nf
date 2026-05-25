#!/usr/bin/env python3
"""
filter_high_contam.py
---------------------
Filter co-assembled clusters with high contamination from cluster_data_updated.json.
Outputs one JSON per passing cluster for downstream GTDB-Tk / final reporting.
"""

import argparse
import json
import logging
import sys
from pathlib import Path


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

DEFAULTS = {
    "max_contamination": 10.0,
    "min_completeness":   0.0,   # 0 = disabled
}


def run(
    updated_json: Path,
    out_dir:      Path,
    out_log:      Path,
    cfg:          dict,
) -> None:

    out_dir.mkdir(parents=True, exist_ok=True)

    with updated_json.open() as fh:
        data = json.load(fh)

    passed, failed = [], []

    for cluster in data["clusters"]:
        cid = cluster["cluster_id"]
        r1  = cluster.get("co_assembly_results", {}).get("round_1", {})

        comp = r1.get("completeness",  0)
        cont = r1.get("contamination", 0)

        fail_reason = None
        if cont > cfg["max_contamination"]:
            fail_reason = f"high_contamination({cont:.2f}%)"
        elif cfg["min_completeness"] > 0 and comp < cfg["min_completeness"]:
            fail_reason = f"low_completeness({comp:.2f}%)"

        cluster["filter_result"] = "fail:" + fail_reason if fail_reason else "pass"

        if fail_reason:
            failed.append(cid)
            log.info("Cluster %s FAIL: %s", cid, fail_reason)
        else:
            passed.append(cid)
            out_path = out_dir / f"cluster_{cid}.json"
            out_path.write_text(json.dumps(cluster, indent=2))
            log.info("Cluster %s PASS: comp=%.1f%% cont=%.1f%%", cid, comp, cont)

    total = len(passed) + len(failed)
    log_lines = [
        f"total_clusters:   {total}",
        f"passed:           {len(passed)}",
        f"failed:           {len(failed)}",
        f"max_contamination:{cfg['max_contamination']}",
        f"min_completeness: {cfg['min_completeness']} (0=disabled)",
        "",
        "# Failed clusters",
    ] + [f"  {cid}" for cid in failed]

    out_log.write_text("\n".join(log_lines) + "\n")

    if not passed:
        log.error("No clusters passed contamination filter.")
        sys.exit(1)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--updated_json",      required=True, type=Path)
    p.add_argument("--out_dir",           required=True, type=Path)
    p.add_argument("--out_log",           required=True, type=Path)
    p.add_argument("--max_contamination", type=float,
                   default=DEFAULTS["max_contamination"])
    p.add_argument("--min_completeness",  type=float,
                   default=DEFAULTS["min_completeness"])
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    cfg  = {k: getattr(args, k) for k in DEFAULTS}
    run(args.updated_json, args.out_dir, args.out_log, cfg)
