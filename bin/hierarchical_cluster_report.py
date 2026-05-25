#!/usr/bin/env python3
"""
hierarchical_cluster_report.py
-------------------------------
Generate final clustering summary and identify high-quality clusters
based on size and internal similarity thresholds.
"""

import argparse
import csv
import time
import numpy as np
from collections import Counter
from pathlib import Path


def load_clusters(path: Path) -> list:
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader)


def load_similarity_csv(path: Path):
    with path.open() as fh:
        reader = csv.reader(fh)
        header = next(reader)
        labels = header[1:]
        matrix = []
        for row in reader:
            matrix.append([float(v) for v in row[1:]])
    label_index = {l: i for i, l in enumerate(labels)}
    return label_index, np.array(matrix)


def run(
    clusters_tsv:            Path,
    similarity_csv:          Path,
    out_summary:             Path,
    out_statistics:          Path,
    out_high_quality:        Path,
    min_cluster_size:        int,
    max_cluster_size:        int,
    min_similarity_in_cluster: float,
) -> None:

    rows = load_clusters(clusters_tsv)
    label_index, sim_matrix = load_similarity_csv(similarity_csv)

    n_sags = len(rows)
    print(f"Loaded {n_sags} SAG cluster assignments")

    # Group SAGs by cluster
    cluster_members: dict = {}
    for row in rows:
        cid = int(row["Cluster_ID"])
        cluster_members.setdefault(cid, []).append(row["SAG_ID"])

    # Per-cluster statistics
    cluster_stats = []
    for cid, members in sorted(cluster_members.items()):
        size = len(members)

        # Internal similarity
        sims = []
        for i, s1 in enumerate(members):
            for j, s2 in enumerate(members):
                if i < j and s1 in label_index and s2 in label_index:
                    sims.append(sim_matrix[label_index[s1], label_index[s2]])

        if sims:
            mean_sim = float(np.mean(sims))
            min_sim  = float(np.min(sims))
            max_sim  = float(np.max(sims))
            std_sim  = float(np.std(sims))
        elif size == 1:
            mean_sim = min_sim = max_sim = std_sim = 1.0   # singleton
        else:
            mean_sim = min_sim = max_sim = std_sim = float("nan")

        # Quality flags
        size_ok = min_cluster_size <= size <= max_cluster_size
        sim_ok  = np.isnan(mean_sim) or mean_sim >= min_similarity_in_cluster
        is_hq   = size_ok and sim_ok

        cluster_stats.append({
            "Cluster_ID":       cid,
            "Size":             size,
            "Mean_Similarity":  mean_sim,
            "Min_Similarity":   min_sim,
            "Max_Similarity":   max_sim,
            "Std_Similarity":   std_sim,
            "High_Quality":     is_hq,
            "Members_Preview":  ",".join(members[:5]) + ("..." if size > 5 else ""),
        })

    # Sort by size descending
    cluster_stats.sort(key=lambda x: -x["Size"])

    # Write cluster statistics TSV
    with out_statistics.open("w", newline="") as fh:
        fieldnames = ["Cluster_ID", "Size", "Mean_Similarity", "Min_Similarity",
                      "Max_Similarity", "Std_Similarity", "High_Quality", "Members_Preview"]
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in cluster_stats:
            writer.writerow({k: (f"{row[k]:.6f}" if isinstance(row[k], float) else row[k])
                             for k in fieldnames})

    # Write high-quality clusters TSV
    hq = [r for r in cluster_stats if r["High_Quality"]]
    with out_high_quality.open("w", newline="") as fh:
        fieldnames = ["Cluster_ID", "Size", "Mean_Similarity",
                      "Min_Similarity", "Max_Similarity", "Members_Preview"]
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in hq:
            writer.writerow({k: (f"{row[k]:.6f}" if isinstance(row[k], float) else row[k])
                             for k in fieldnames})
    print(f"High-quality clusters: {len(hq)}/{len(cluster_stats)}")

    # Summary
    n_clusters  = len(cluster_stats)
    sizes       = [r["Size"] for r in cluster_stats]
    size_dist   = Counter(sizes)
    singletons  = size_dist.get(1, 0)
    multi       = sum(v for k, v in size_dist.items() if k > 1)
    large_warn  = sum(1 for s in sizes if s > max_cluster_size)

    with out_summary.open("w") as fh:
        fh.write("# Final Hierarchical Clustering Summary\n")
        fh.write(f"generated_on:        {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        fh.write("## Overall Statistics\n")
        fh.write(f"total_sags:          {n_sags}\n")
        fh.write(f"total_clusters:      {n_clusters}\n")
        fh.write(f"singleton_clusters:  {singletons}\n")
        fh.write(f"multi_sag_clusters:  {multi}\n")
        fh.write(f"high_quality:        {len(hq)}\n")
        fh.write(f"large_cluster_warn:  {large_warn} (>{max_cluster_size} SAGs)\n\n")
        fh.write("## Filters Applied\n")
        fh.write(f"min_cluster_size:           {min_cluster_size}\n")
        fh.write(f"max_cluster_size:           {max_cluster_size}\n")
        fh.write(f"min_similarity_in_cluster:  {min_similarity_in_cluster}\n\n")
        fh.write("## Cluster Size Distribution\n")
        for sz in sorted(size_dist.keys()):
            fh.write(f"size_{sz:03d}: {size_dist[sz]}\n")
        fh.write("\n## Top 10 Largest Clusters\n")
        for r in cluster_stats[:10]:
            fh.write(f"Cluster {r['Cluster_ID']:4d}: "
                     f"{r['Size']:3d} SAGs, "
                     f"mean_sim={r['Mean_Similarity']:.4f}\n")
        fh.write("\n## High-Quality Clusters\n")
        for r in hq:
            fh.write(f"Cluster {r['Cluster_ID']:4d}: "
                     f"{r['Size']:3d} SAGs, "
                     f"mean_sim={r['Mean_Similarity']:.4f} "
                     f"± {r['Std_Similarity']:.4f}\n")
        fh.write("\n## Recommendations\n")
        fh.write("1. Use high_quality_clusters.tsv for co-assembly input\n")
        fh.write("2. Singleton clusters represent unique genomes; skip co-assembly\n")
        if large_warn:
            fh.write(f"3. WARNING: {large_warn} clusters exceed max_cluster_size={max_cluster_size}; "
                     "consider stricter threshold\n")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--clusters_tsv",              required=True, type=Path)
    p.add_argument("--similarity_csv",            required=True, type=Path)
    p.add_argument("--out_summary",               required=True, type=Path)
    p.add_argument("--out_statistics",            required=True, type=Path)
    p.add_argument("--out_high_quality",          required=True, type=Path)
    p.add_argument("--min_cluster_size",          type=int,   default=2)
    p.add_argument("--max_cluster_size",          type=int,   default=50)
    p.add_argument("--min_similarity_in_cluster", type=float, default=0.2)
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(
        clusters_tsv              = args.clusters_tsv,
        similarity_csv            = args.similarity_csv,
        out_summary               = args.out_summary,
        out_statistics            = args.out_statistics,
        out_high_quality          = args.out_high_quality,
        min_cluster_size          = args.min_cluster_size,
        max_cluster_size          = args.max_cluster_size,
        min_similarity_in_cluster = args.min_similarity_in_cluster,
    )