#!/usr/bin/env python3
"""
hierarchical_cluster_analysis.py
---------------------------------
Perform hierarchical clustering on a distance matrix and generate
cluster assignments, dendrogram, and validation metrics.
"""

import argparse
import csv
import sys
import time
import numpy as np
from collections import Counter
from pathlib import Path


def load_distance_csv(path: Path):
    with path.open() as fh:
        reader = csv.reader(fh)
        header = next(reader)
        labels = header[1:]
        matrix = []
        for row in reader:
            matrix.append([float(v) for v in row[1:]])
    return labels, np.array(matrix)


def run(
    distance_csv:        Path,
    out_clusters:        Path,
    out_linkage:         Path,
    out_dendrogram:      Path,
    out_report:          Path,
    out_validation:      Path,
    linkage_method:      str,
    criterion:           str,
    threshold:           float,
    generate_dendrogram: bool,
    compute_validation:  bool,
    min_cluster_size:    int,
    max_cluster_size:    int,
) -> None:

    import scipy
    from scipy.spatial.distance import squareform
    from scipy.cluster.hierarchy import (
        linkage, fcluster, dendrogram as scipy_dendrogram,
        inconsistent, cophenet
    )

    labels, dist_matrix = load_distance_csv(distance_csv)
    n = len(labels)
    print(f"Loaded {n}x{n} distance matrix")
    print(f"Linkage: {linkage_method}, Criterion: {criterion}, Threshold: {threshold}")

    # Ensure symmetry and non-negative
    dist_matrix = (dist_matrix + dist_matrix.T) / 2
    np.fill_diagonal(dist_matrix, 0.0)
    dist_matrix = np.clip(dist_matrix, 0.0, None)

    # Condensed distance vector
    condensed = squareform(dist_matrix, checks=False)

    # Hierarchical clustering
    t0 = time.time()
    Z = linkage(condensed, method=linkage_method)
    elapsed = time.time() - t0
    print(f"Linkage computed in {elapsed:.2f}s")

    # Form flat clusters
    if criterion == "maxclust":
        clusters = fcluster(Z, int(threshold), criterion=criterion)
    else:
        clusters = fcluster(Z, threshold, criterion=criterion)

    cluster_counts = Counter(clusters)
    n_clusters = len(cluster_counts)
    print(f"Clusters formed: {n_clusters}")

    # Write cluster assignments
    with out_clusters.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["SAG_ID", "Cluster_ID", "Cluster_Size"])
        for sag, cid in sorted(zip(labels, clusters), key=lambda x: (x[1], x[0])):
            writer.writerow([sag, cid, cluster_counts[cid]])
    print(f"Cluster assignments written: {len(labels)} SAGs")

    # Write linkage matrix
    with out_linkage.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["SAG1_Index", "SAG2_Index", "Distance", "Cluster_Size"])
        for row in Z:
            writer.writerow([f"{v:.6f}" for v in row])

    # Dendrogram
    if generate_dendrogram:
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt

            fig_width = min(max(12, n * 0.3), 200)
            plt.figure(figsize=(fig_width, 8))
            scipy_dendrogram(
                Z,
                labels=labels if n <= 50 else None,
                leaf_rotation=90,
                leaf_font_size=6 if n > 50 else 8,
            )
            plt.title(
                f"Hierarchical Clustering Dendrogram ({n} SAGs)\n"
                f"Method: {linkage_method}, Criterion: {criterion}, "
                f"Threshold: {threshold}"
            )
            plt.xlabel("SAG ID")
            plt.ylabel("Distance")
            plt.tight_layout()
            plt.savefig(str(out_dendrogram), dpi=300, bbox_inches='tight')
            plt.close()
            print(f"Dendrogram saved: {out_dendrogram}")
        except Exception as e:
            print(f"WARNING: Could not generate dendrogram: {e}")
            out_dendrogram.write_text("# Dendrogram generation failed\n")
    else:
        out_dendrogram.write_text("# Dendrogram generation disabled\n")

    # Validation metrics
    cophenetic_corr = np.nan
    mean_within     = np.nan
    mean_between    = np.nan

    if compute_validation:
        try:
            c, _ = cophenet(Z, condensed)
            cophenetic_corr = float(c)
            print(f"Cophenetic correlation: {cophenetic_corr:.4f}")
        except Exception as e:
            print(f"WARNING: Cophenetic correlation failed: {e}")

        within, between = [], []
        for i in range(n):
            for j in range(i + 1, n):
                d = dist_matrix[i, j]
                if clusters[i] == clusters[j]:
                    within.append(d)
                else:
                    between.append(d)
        mean_within  = float(np.mean(within))  if within  else 0.0
        mean_between = float(np.mean(between)) if between else 0.0
        print(f"Mean within-cluster distance:  {mean_within:.4f}")
        print(f"Mean between-cluster distance: {mean_between:.4f}")

    # Cluster quality flags
    large_clusters   = [cid for cid, sz in cluster_counts.items() if sz > max_cluster_size]
    small_clusters   = [cid for cid, sz in cluster_counts.items() if sz < min_cluster_size]
    if large_clusters:
        print(f"WARNING: {len(large_clusters)} clusters exceed max_cluster_size={max_cluster_size}")

    # Clustering report
    size_dist = Counter(cluster_counts.values())
    with out_report.open("w") as fh:
        fh.write("# Hierarchical Clustering Analysis Report\n")
        fh.write(f"scipy_version:     {scipy.__version__}\n")
        fh.write(f"linkage_method:    {linkage_method}\n")
        fh.write(f"criterion:         {criterion}\n")
        fh.write(f"threshold:         {threshold}\n")
        fh.write(f"n_sags:            {n}\n")
        fh.write(f"n_clusters:        {n_clusters}\n")
        fh.write(f"clustering_time_s: {elapsed:.2f}\n\n")
        fh.write("## Cluster Size Distribution\n")
        for size in sorted(size_dist.keys()):
            fh.write(f"size_{size:03d}: {size_dist[size]} clusters\n")
        fh.write(f"\n## Quality Flags\n")
        fh.write(f"clusters_exceeding_max_size ({max_cluster_size}): {len(large_clusters)}\n")
        fh.write(f"clusters_below_min_size ({min_cluster_size}):     {len(small_clusters)}\n")
        fh.write("\n## Detailed Cluster Assignments\n")
        for cid in sorted(cluster_counts.keys()):
            members = [s for s, c in zip(labels, clusters) if c == cid]
            preview = ", ".join(members[:10])
            extra   = f" ... +{len(members)-10} more" if len(members) > 10 else ""
            fh.write(f"Cluster {cid:4d} ({len(members):3d} SAGs): {preview}{extra}\n")

    # Validation report
    with out_validation.open("w") as fh:
        fh.write("# Cluster Validation Report\n")
        if np.isnan(cophenetic_corr):
            fh.write("cophenetic_correlation: NA\n")
            fh.write("interpretation:         not computed\n")
        else:
            quality = (
                "excellent" if cophenetic_corr > 0.9 else
                "good"      if cophenetic_corr > 0.8 else
                "fair"      if cophenetic_corr > 0.7 else
                "poor"
            )
            fh.write(f"cophenetic_correlation: {cophenetic_corr:.6f}\n")
            fh.write(f"interpretation:         {quality}\n")
        fh.write(f"\nmean_within_cluster_distance:  {mean_within:.6f}\n")
        fh.write(f"mean_between_cluster_distance: {mean_between:.6f}\n")
        if mean_within and mean_within > 0:
            fh.write(f"separation_ratio:              {mean_between/mean_within:.6f}\n")
        else:
            fh.write("separation_ratio:              NA\n")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--distance_csv",        required=True,  type=Path)
    p.add_argument("--out_clusters",        required=True,  type=Path)
    p.add_argument("--out_linkage",         required=True,  type=Path)
    p.add_argument("--out_dendrogram",      required=True,  type=Path)
    p.add_argument("--out_report",          required=True,  type=Path)
    p.add_argument("--out_validation",      required=True,  type=Path)
    p.add_argument("--linkage_method",      default="complete",
                   choices=["complete", "average", "single", "ward"])
    p.add_argument("--criterion",           default="inconsistent",
                   choices=["inconsistent", "distance", "maxclust"])
    p.add_argument("--threshold",           type=float, default=0.95)
    p.add_argument("--generate_dendrogram", action="store_true", default=True)
    p.add_argument("--no_dendrogram",       dest="generate_dendrogram",
                   action="store_false")
    p.add_argument("--compute_validation",  action="store_true", default=True)
    p.add_argument("--no_validation",       dest="compute_validation",
                   action="store_false")
    p.add_argument("--min_cluster_size",    type=int, default=2)
    p.add_argument("--max_cluster_size",    type=int, default=50)
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(
        distance_csv        = args.distance_csv,
        out_clusters        = args.out_clusters,
        out_linkage         = args.out_linkage,
        out_dendrogram      = args.out_dendrogram,
        out_report          = args.out_report,
        out_validation      = args.out_validation,
        linkage_method      = args.linkage_method,
        criterion           = args.criterion,
        threshold           = args.threshold,
        generate_dendrogram = args.generate_dendrogram,
        compute_validation  = args.compute_validation,
        min_cluster_size    = args.min_cluster_size,
        max_cluster_size    = args.max_cluster_size,
    )
