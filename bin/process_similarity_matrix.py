#!/usr/bin/env python3
"""
process_similarity_matrix.py
-----------------------------
Generate statistics and TSV report from sourmash compare output.
Note: sourmash compare --csv already produces the labeled CSV directly.
This script adds statistics reporting on top.
"""

import argparse
import csv
import time
import numpy as np
from pathlib import Path


def load_labels(labels_file: Path) -> list:
    """
    Load SAG IDs from sourmash-generated labels file.
    Labels are exactly the --name values used during sketch,
    which are already SAG IDs (e.g. AACAACACGAAACCACCACGH).
    """
    labels = []
    with labels_file.open() as fh:
        for line in fh:
            name = line.strip()
            if name:
                labels.append(name)
    return labels


def run(
    matrix_npy:  Path,
    labels_file: Path,
    matrix_csv:  Path,   # already produced by sourmash compare --csv
    ksize:       int,
    scaled:      int,
    out_report:  Path,
) -> None:

    t0 = time.time()

    # Load matrix and labels
    matrix = np.load(str(matrix_npy))
    labels = load_labels(labels_file)
    n = len(labels)

    assert matrix.shape == (n, n), (
        f"Matrix shape {matrix.shape} does not match {n} labels"
    )

    print(f"Loaded {n}x{n} similarity matrix")
    print(f"Sample labels: {labels[:3]}")

    # Statistics on off-diagonal values
    off_diag = matrix[~np.eye(n, dtype=bool)]
    n_pairs  = n * (n - 1) // 2

    stats = {
        "mean":   float(np.mean(off_diag)),
        "median": float(np.median(off_diag)),
        "max":    float(np.max(off_diag)),
        "min":    float(np.min(off_diag)),
        "std":    float(np.std(off_diag)),
    }

    # Count pairs above common thresholds
    thresholds = [0.01, 0.05, 0.1, 0.2, 0.5]
    threshold_counts = {
        t: int(np.sum(off_diag >= t)) // 2
        for t in thresholds
    }

    elapsed = time.time() - t0

    with out_report.open("w") as fh:
        fh.write("# Sourmash Similarity Matrix Report\n")
        fh.write(f"ksize:             {ksize}\n")
        fh.write(f"scaled:            {scaled}\n")
        fh.write(f"n_sags:            {n}\n")
        fh.write(f"n_pairs:           {n_pairs}\n")
        fh.write(f"mean_similarity:   {stats['mean']:.6f}\n")
        fh.write(f"median_similarity: {stats['median']:.6f}\n")
        fh.write(f"max_similarity:    {stats['max']:.6f}\n")
        fh.write(f"min_similarity:    {stats['min']:.6f}\n")
        fh.write(f"std_similarity:    {stats['std']:.6f}\n")
        fh.write(f"elapsed_sec:       {elapsed:.1f}\n")
        fh.write("\n# Pairs above similarity thresholds\n")
        for t, count in threshold_counts.items():
            fh.write(f"pairs_above_{t}: {count}\n")

    print(f"Report written.")
    print(f"Mean similarity: {stats['mean']:.6f}")
    print(f"Max similarity:  {stats['max']:.6f}")
    for t, count in threshold_counts.items():
        print(f"  Pairs >= {t}: {count}")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--matrix_npy",  required=True, type=Path)
    p.add_argument("--labels",      required=True, type=Path)
    p.add_argument("--matrix_csv",  required=True, type=Path,
                   help="CSV already produced by sourmash compare --csv")
    p.add_argument("--ksize",       required=True, type=int)
    p.add_argument("--scaled",      required=True, type=int)
    p.add_argument("--out_report",  required=True, type=Path)
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(
        matrix_npy  = args.matrix_npy,
        labels_file = args.labels,
        matrix_csv  = args.matrix_csv,
        ksize       = args.ksize,
        scaled      = args.scaled,
        out_report  = args.out_report,
    )