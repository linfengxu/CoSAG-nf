#!/usr/bin/env python3
"""
Build Round 1 CoSAG JSON and merge reads per cluster.

Usage:
    prepare_cluster_json.py \
        --samples   paired_end_qc_samples.tsv \
        --clusters  hierarchical_clusters.tsv \
        --contigs   sag_contigs_mapping.tsv \
        --checkm2   sags_pass.tsv \
        --outdir    . \
        --out       round1_clusters.json \
        --threads   8
"""

import argparse
import json
import logging
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ─────────────────────────────────────────────
# Parsers
# ─────────────────────────────────────────────

def parse_samples(path):
    result = {}
    with open(path) as f:
        f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            result[parts[0]] = {"read1": parts[1], "read2": parts[2]}
    log.info("  samples:  %d SAGs", len(result))
    return result


def parse_clusters(path):
    result = defaultdict(list)
    with open(path) as f:
        f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            result[parts[1]].append(parts[0])
    log.info("  clusters: %d clusters", len(result))
    return dict(result)


def parse_contigs(path):
    result = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            result[parts[0]] = parts[1]
    log.info("  contigs:  %d SAGs", len(result))
    return result


def parse_checkm2(path):
    result = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {col: i for i, col in enumerate(header)}
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if not parts or len(parts) < len(idx):
                continue
            name = parts[idx["Name"]]
            sag_id = name[:-8] if name.endswith("_contigs") else name
            comp = float(parts[idx["Completeness"]])
            cont = float(parts[idx["Contamination"]])
            result[sag_id] = {
                "completeness":   comp,
                "contamination":  cont,
                "quality_score":  round(comp - 5 * cont, 3),
                "genome_size":    int(parts[idx["Genome_Size"]]),
                "contig_n50":     int(parts[idx["Contig_N50"]]),
                "gc_content":     float(parts[idx["GC_Content"]]),
                "total_contigs":  int(parts[idx["Total_Contigs"]]),
                "coding_density": float(parts[idx["Coding_Density"]]),
            }
    log.info("  checkm2:  %d SAGs", len(result))
    return result


# ─────────────────────────────────────────────
# Read merging
# ─────────────────────────────────────────────

def merge_one_cluster(cluster_id, r1_list, r2_list, outdir):
    """Merge R1 and R2 reads for a single cluster. Returns (cluster_id, ok, err)."""
    cluster_dir = (Path(outdir) / "cluster_{}".format(cluster_id)).resolve()
    cluster_dir.mkdir(parents=True, exist_ok=True)

    r1_out = cluster_dir / "merged_R1.fastq.gz"
    r2_out = cluster_dir / "merged_R2.fastq.gz"

    for inputs, out in [(r1_list, r1_out), (r2_list, r2_out)]:
        # Filter out None paths
        valid = [p for p in inputs if p]
        if not valid:
            return cluster_id, False, "no valid read paths"
        cmd = "cat {} > {}".format(" ".join('"{}"'.format(p) for p in valid), out)
        ret = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        if ret.returncode != 0:
            return cluster_id, False, ret.stderr.decode().strip()

    return cluster_id, True, None


def merge_all_clusters(cluster_reads, outdir, threads):
    """
    cluster_reads: {cluster_id: {"r1": [...], "r2": [...]}}
    Returns: {cluster_id: bool}  True = success
    """
    results = {}
    log.info("Merging reads for %d clusters (threads=%d)...", len(cluster_reads), threads)

    with ThreadPoolExecutor(max_workers=threads) as pool:
        futures = {
            pool.submit(merge_one_cluster, cid, data["r1"], data["r2"], outdir): cid
            for cid, data in cluster_reads.items()
        }
        for future in as_completed(futures):
            cid, ok, err = future.result()
            results[cid] = ok
            if ok:
                log.info("  [ok] cluster_%s", cid)
            else:
                log.error("  [fail] cluster_%s: %s", cid, err)

    n_ok = sum(results.values())
    log.info("Merge done: %d/%d succeeded", n_ok, len(results))
    if n_ok < len(results):
        log.error("%d clusters failed to merge", len(results) - n_ok)
        sys.exit(1)

    return results


# ─────────────────────────────────────────────
# JSON builder
# ─────────────────────────────────────────────

def build_json(samples, clusters, contigs, checkm2, outdir):
    missing_reads   = []
    missing_contigs = []
    missing_checkm2 = []

    cluster_list  = []
    cluster_reads = {}   # for merge step

    for cluster_id in sorted(clusters.keys(), key=lambda x: int(x) if x.isdigit() else x):
        sag_ids = clusters[cluster_id]
        members = []
        r1_list = []
        r2_list = []

        for sag_id in sag_ids:
            member = {"sag_id": sag_id}

            if sag_id in samples:
                member["read1"] = samples[sag_id]["read1"]
                member["read2"] = samples[sag_id]["read2"]
                r1_list.append(samples[sag_id]["read1"])
                r2_list.append(samples[sag_id]["read2"])
            else:
                member["read1"] = None
                member["read2"] = None
                missing_reads.append(sag_id)

            if sag_id in contigs:
                member["individual_contigs"] = contigs[sag_id]
            else:
                member["individual_contigs"] = None
                missing_contigs.append(sag_id)

            if sag_id in checkm2:
                member["checkm2"] = checkm2[sag_id]
            else:
                member["checkm2"] = None
                missing_checkm2.append(sag_id)

            members.append(member)

        cluster_dir = (Path(outdir) / "cluster_{}".format(cluster_id)).resolve()
        cluster_reads[cluster_id] = {"r1": r1_list, "r2": r2_list}

        cluster_list.append({
            "cluster_id":         "R1_C{}".format(cluster_id),
            "cluster_size":       len(sag_ids),
            "leaf_sags":          sag_ids,
            "members":            members,
            "merged_reads": {
                "read1": str(cluster_dir / "merged_R1.fastq.gz"),
                "read2": str(cluster_dir / "merged_R2.fastq.gz"),
            },
            "coassembly_contigs": None,
            "checkm2_coassembly": None,
            "status":             "pending",
        })

    if missing_reads:
        log.warning("  %d SAGs missing reads (first 5): %s", len(missing_reads), missing_reads[:5])
    if missing_contigs:
        log.warning("  %d SAGs missing contigs (first 5): %s", len(missing_contigs), missing_contigs[:5])
    if missing_checkm2:
        log.warning("  %d SAGs missing checkm2 (first 5): %s", len(missing_checkm2), missing_checkm2[:5])

    data = {
        "version":     "1.0",
        "round":       1,
        "member_type": "sag",
        "n_clusters":  len(cluster_list),
        "n_sags":      sum(len(clusters[c]) for c in clusters),
        "clusters":    cluster_list,
    }

    return data, cluster_reads


# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples",  required=True)
    parser.add_argument("--clusters", required=True)
    parser.add_argument("--contigs",  required=True)
    parser.add_argument("--checkm2",  required=True)
    parser.add_argument("--outdir",   required=True)
    parser.add_argument("--out",      default="round1_clusters.json")
    parser.add_argument("--threads",  type=int, default=8)
    args = parser.parse_args()

    log.info("Loading input files...")
    samples  = parse_samples(args.samples)
    clusters = parse_clusters(args.clusters)
    contigs  = parse_contigs(args.contigs)
    checkm2  = parse_checkm2(args.checkm2)

    log.info("Building JSON structure...")
    data, cluster_reads = build_json(samples, clusters, contigs, checkm2, args.outdir)

    log.info("Merging reads...")
    merge_all_clusters(cluster_reads, args.outdir, args.threads)

    out_path = Path(args.out)
    with open(str(out_path), "w") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)

    log.info("Done -> %s  (clusters: %d, SAGs: %d)", out_path, data["n_clusters"], data["n_sags"])


if __name__ == "__main__":
    main()
