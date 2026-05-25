#!/usr/bin/env python3
"""
parse_cluster_json.py
---------------------
Parse cluster_data.json and emit one-line-per-cluster TSV
for Nextflow channel creation.

Output TSV columns:
    cluster_id  read1_merged_path  read2_merged_path  member_count  member_ids
"""

import argparse
import json
import sys
from pathlib import Path


def run(cluster_json: Path, out_tsv: Path, status_filter: list) -> None:

    with cluster_json.open() as fh:
        data = json.load(fh)

    rows = []
    for cluster in data["clusters"]:
        cid   = cluster["cluster_id"]
        hq    = cluster.get("high_quality", True)

        # collect pass members only
        pass_members = [
            m for m in cluster["members"]
            if m.get("status") in status_filter
        ]

        if not pass_members:
            print(f"Cluster {cid}: no passing members, skipping", file=sys.stderr)
            continue

        if not hq:
            print(f"Cluster {cid}: not high-quality, skipping", file=sys.stderr)
            continue

        read1s  = [m["read1"]  for m in pass_members if m.get("read1")]
        read2s  = [m["read2"]  for m in pass_members if m.get("read2")]
        members = [m["sag_id"] for m in pass_members]

        rows.append("\t".join([
            cid,
            ",".join(read1s),
            ",".join(read2s),
            str(len(members)),
            ",".join(members),
        ]))

    out_tsv.write_text("\n".join(rows) + "\n")
    print(f"Wrote {len(rows)} clusters to {out_tsv}")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cluster_json", required=True, type=Path)
    p.add_argument("--out_tsv",      required=True, type=Path)
    p.add_argument("--status_filter", nargs="+",
                   default=["pass"],
                   help="SAG statuses to include (default: pass)")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(args.cluster_json, args.out_tsv, args.status_filter)
