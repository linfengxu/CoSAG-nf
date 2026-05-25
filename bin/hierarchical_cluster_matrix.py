#!/usr/bin/env python3
"""
hierarchical_cluster_matrix.py
-------------------------------
Convert Jaccard similarity matrix to distance matrix for hierarchical clustering.

Conversion formulas:
    euclidean : distance = sqrt(2 * (1 - similarity))
    cosine    : distance = 1 - similarity
    manhattan : distance = 2 * (1 - similarity)
"""

import argparse
import csv
import sys
import numpy as np
from pathlib import Path


DISTANCE_FORMULAS = {
    "euclidean": lambda s: np.sqrt(2 * (1 - s)),
    "cosine":    lambda s: 1 - s,
    "manhattan": lambda s: 2 * (1 - s),
}


def load_similarity_csv(path: Path):
    """Load sourmash compare --csv output."""
    with path.open() as fh:
        reader = csv.reader(fh)
        header = next(reader)
        col_labels = header[1:]          # 从表头行提取列标签

        row_labels = []
        matrix = []
        for row in reader:
            if not row or all(v.strip() == "" for v in row):
                continue                  # 跳过末尾空行
            row_labels.append(row[0])
            matrix.append([float(v) for v in row[1:]])

    # 检查行列标签是否一致
    if col_labels != row_labels:
        print(
            f"WARNING: 列标签数 ({len(col_labels)}) 与行标签数 "
            f"({len(row_labels)}) 不一致，取最小值并裁剪矩阵。"
        )

    n = min(len(col_labels), len(row_labels), len(matrix))
    labels = row_labels[:n]
    mat = np.array([r[:n] for r in matrix[:n]])
    return labels, mat


def run(
    similarity_csv: Path,
    out_distance:   Path,
    out_stats:      Path,
    distance_metric: str,
) -> None:

    if not similarity_csv.exists():
        print(f"ERROR: Similarity matrix not found: {similarity_csv}")
        sys.exit(1)

    labels, sim_matrix = load_similarity_csv(similarity_csv)
    n = len(labels)
    print(f"Loaded {n}x{n} similarity matrix")

    if n < 2:
        print(f"ERROR: Need at least 2 SAGs for clustering, got {n}")
        sys.exit(1)

    # Convert similarity to distance
    if distance_metric not in DISTANCE_FORMULAS:
        print(f"WARNING: Unknown metric '{distance_metric}', using euclidean")
        distance_metric = "euclidean"

    dist_matrix = DISTANCE_FORMULAS[distance_metric](sim_matrix)

    # Fix numerical issues
    dist_matrix = np.nan_to_num(dist_matrix, nan=0.0, posinf=np.nanmax(dist_matrix[np.isfinite(dist_matrix)]))
    np.fill_diagonal(dist_matrix, 0.0)
    dist_matrix = (dist_matrix + dist_matrix.T) / 2   # enforce symmetry
    dist_matrix = np.clip(dist_matrix, 0.0, None)      # enforce non-negative

    # Write distance matrix CSV
    with out_distance.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([""] + labels)
        for i, name in enumerate(labels):
            writer.writerow([name] + [f"{v:.6f}" for v in dist_matrix[i]])
    print(f"Distance matrix written: {n}x{n}")

    # Statistics
    off_diag_sim  = sim_matrix[~np.eye(n, dtype=bool)]
    off_diag_dist = dist_matrix[~np.eye(n, dtype=bool)]

    with out_stats.open("w") as fh:
        fh.write("# Similarity → Distance Conversion Report\n")
        fh.write(f"distance_metric:    {distance_metric}\n")
        fh.write(f"n_sags:             {n}\n\n")
        fh.write("## Similarity Statistics\n")
        fh.write(f"mean:   {np.mean(off_diag_sim):.6f}\n")
        fh.write(f"median: {np.median(off_diag_sim):.6f}\n")
        fh.write(f"std:    {np.std(off_diag_sim):.6f}\n")
        fh.write(f"min:    {np.min(off_diag_sim):.6f}\n")
        fh.write(f"max:    {np.max(off_diag_sim):.6f}\n\n")
        fh.write("## Distance Statistics\n")
        fh.write(f"mean:   {np.mean(off_diag_dist):.6f}\n")
        fh.write(f"median: {np.median(off_diag_dist):.6f}\n")
        fh.write(f"std:    {np.std(off_diag_dist):.6f}\n")
        fh.write(f"min:    {np.min(off_diag_dist):.6f}\n")
        fh.write(f"max:    {np.max(off_diag_dist):.6f}\n\n")
        fh.write("## Matrix Properties\n")
        fh.write(f"symmetric:      {np.allclose(dist_matrix, dist_matrix.T)}\n")
        fh.write(f"zero_diagonal:  {np.allclose(np.diag(dist_matrix), 0)}\n")
        fh.write(f"non_negative:   {bool(np.all(dist_matrix >= 0))}\n")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--similarity_csv",  required=True, type=Path)
    p.add_argument("--out_distance",    required=True, type=Path)
    p.add_argument("--out_stats",       required=True, type=Path)
    p.add_argument("--distance_metric", default="euclidean",
                   choices=list(DISTANCE_FORMULAS.keys()))
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(
        similarity_csv  = args.similarity_csv,
        out_distance    = args.out_distance,
        out_stats       = args.out_stats,
        distance_metric = args.distance_metric,
    )
