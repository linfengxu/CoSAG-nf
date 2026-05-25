#!/usr/bin/env python3
import argparse
import json
from copy import deepcopy
from pathlib import Path


def load_json(path: Path):
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def build_round1_lookup(round1_data: dict):
    lookup = {}
    for cluster in round1_data.get("clusters", []):
        cid = cluster.get("cluster_id")
        if not cid:
            continue
        lookup[cid] = {
            "gtdbtk_classification": deepcopy(cluster.get("co_assembly_results", {}).get("gtdbtk_classification")),
            "leaf_gtdbtk": deepcopy(cluster.get("leaf_gtdbtk")),
        }
    return lookup


def backfill(round2_data: dict, lookup: dict):
    touched = 0
    for cluster in round2_data.get("clusters", []):
        for member in cluster.get("members", []):
            cid = member.get("cluster_id")
            if not cid or cid not in lookup:
                continue
            src = lookup[cid]
            if member.get("gtdbtk_classification") in (None, {}, "") and src.get("gtdbtk_classification"):
                member["gtdbtk_classification"] = deepcopy(src["gtdbtk_classification"])
                touched += 1
            if (not member.get("leaf_gtdbtk")) and src.get("leaf_gtdbtk"):
                member["leaf_gtdbtk"] = deepcopy(src["leaf_gtdbtk"])
                touched += 1
    return touched


def main():
    parser = argparse.ArgumentParser(
        description="Backfill round2 member GTDB annotations from round1 final JSON."
    )
    parser.add_argument("--round1-final-json", required=True, help="Path to round1 cluster_data_gtdbtk.json")
    parser.add_argument("--round2-updated-json", required=True, help="Path to round2 updated_clusters.json")
    parser.add_argument("--out-json", default=None, help="Output JSON path (default: overwrite round2-updated-json)")
    args = parser.parse_args()

    round1_path = Path(args.round1_final_json)
    round2_path = Path(args.round2_updated_json)
    out_path = Path(args.out_json) if args.out_json else round2_path

    round1_data = load_json(round1_path)
    round2_data = load_json(round2_path)

    lookup = build_round1_lookup(round1_data)
    touched = backfill(round2_data, lookup)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as fh:
        json.dump(round2_data, fh, ensure_ascii=False, indent=2)

    print(f"Backfill complete. Updated fields: {touched}")
    print(f"Wrote: {out_path}")


if __name__ == "__main__":
    main()
