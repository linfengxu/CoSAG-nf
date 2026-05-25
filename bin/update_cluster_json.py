#!/usr/bin/env python3
"""
Update round1_clusters.json with coassembly results.

Inputs:
  --json         round1_clusters.json  (from PREPARE_CLUSTER_JSON)
  --checkm2      quality_report.tsv    (from CHECKM2_COASSEMBLY, batch run)
  --contig-list  coassembly_contig_list.txt  (cluster_id <TAB> /abs/path/contigs.fasta)

Outputs:
  --out          updated_clusters.json         (full updated JSON)
  --per-cluster-dir  ./per_cluster_json/       (one JSON per cluster, for branching)
"""

import argparse
import json
import logging
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


def parse_checkm2(path):
    """
    cluster_id → {completeness, contamination, quality_score, ...}
    Name 列格式: R1_C1_contigs → strip '_contigs' 后得到 cluster_id
    """
    result = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {col: i for i, col in enumerate(header)}
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if not parts or len(parts) < len(idx):
                continue
            name = parts[idx["Name"]]
            # R1_C1_contigs → R1_C1
            cluster_id = name[:-8] if name.endswith("_contigs") else name
            comp  = float(parts[idx["Completeness"]])
            cont  = float(parts[idx["Contamination"]])
            result[cluster_id] = {
                "completeness":   comp,
                "contamination":  cont,
                "quality_score":  round(comp - 5 * cont, 3),
                "genome_size":    int(parts[idx["Genome_Size"]]),
                "contig_n50":     int(parts[idx["Contig_N50"]]),
                "gc_content":     float(parts[idx["GC_Content"]]),
                "total_contigs":  int(parts[idx["Total_Contigs"]]),
                "coding_density": float(parts[idx["Coding_Density"]]),
            }
    log.info("  checkm2: %d clusters", len(result))
    return result


def parse_contig_list(path):
    """cluster_id → /abs/path/contigs.fasta"""
    result = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            result[parts[0]] = parts[1]
    log.info("  contigs: %d clusters", len(result))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--json",            required=True)
    parser.add_argument("--checkm2",         required=True)
    parser.add_argument("--contig-list",     required=True, dest="contig_list")
    parser.add_argument("--out",             default="updated_clusters.json")
    parser.add_argument("--per-cluster-dir", default="per_cluster_json", dest="per_cluster_dir")
    args = parser.parse_args()

    # ── Load inputs ──
    with open(args.json) as f:
        data = json.load(f)

    checkm2     = parse_checkm2(args.checkm2)
    contig_list = parse_contig_list(args.contig_list)

    # ── Update each cluster ──
    missing_checkm2 = []
    missing_contigs = []

    for cluster in data["clusters"]:
        cid = cluster["cluster_id"]

        if cid in contig_list:
            cluster["coassembly_contigs"] = contig_list[cid]
        else:
            missing_contigs.append(cid)
            log.warning("  no contig path for %s", cid)

        if cid in checkm2:
            cluster["checkm2_coassembly"] = checkm2[cid]
            cluster["status"] = "completed"
        else:
            missing_checkm2.append(cid)
            log.warning("  no checkm2 result for %s", cid)

    if missing_contigs:
        log.warning("  %d clusters missing contigs (first 5): %s",
                    len(missing_contigs), missing_contigs[:5])
    if missing_checkm2:
        log.warning("  %d clusters missing checkm2 (first 5): %s",
                    len(missing_checkm2), missing_checkm2[:5])

    # ── Write full updated JSON ──
    out_path = Path(args.out)
    with open(str(out_path), "w") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    log.info("Full JSON → %s", out_path)

    # ── Write per-cluster JSON files ──
    per_dir = Path(args.per_cluster_dir)
    per_dir.mkdir(parents=True, exist_ok=True)

    for cluster in data["clusters"]:
        cid = cluster["cluster_id"]
        out = per_dir / "{}.json".format(cid)
        with open(str(out), "w") as f:
            json.dump(cluster, f, indent=2, ensure_ascii=False)

    log.info("Per-cluster JSON → %s/ (%d files)", per_dir, len(data["clusters"]))


if __name__ == "__main__":
    main()