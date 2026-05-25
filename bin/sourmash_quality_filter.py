#!/usr/bin/env python3
"""
sourmash_quality_filter.py
--------------------------
Optional filter: remove SAGs with insufficient MinHash connections.
Only runs when sourmash_enable_quality_filter = true.
Input: similarity_matrix.csv from sourmash compare --csv
"""

import argparse
import csv
import numpy as np
from pathlib import Path


DEFAULTS = {
    "min_similarity":  0.05,
    "min_connections": 1,
}


def load_matrix_csv(path: Path):
    """Load sourmash compare --csv output. First row and column are labels."""
    with path.open() as fh:
        reader = csv.reader(fh)
        header = next(reader)
        labels = header[1:]     # skip empty first cell
        matrix = []
        for row in reader:
            matrix.append([float(v) for v in row[1:]])
    return labels, np.array(matrix)


def run(
    matrix_csv:     Path,
    out_matrix:     Path,
    out_report:     Path,
    out_low_quality: Path,
    cfg:            dict,
) -> None:

    labels, matrix = load_matrix_csv(matrix_csv)
    n = len(labels)
    print(f"Loaded {n}x{n} matrix")

    # Connectivity analysis
    connectivity = {}
    for i, sag in enumerate(labels):
        row = matrix[i].copy()
        row[i] = 0.0    # exclude self-similarity
        connectivity[sag] = int(np.sum(row >= cfg["min_similarity"]))

    low_quality  = [s for s, c in connectivity.items()
                    if c < cfg["min_connections"]]
    high_quality = [s for s in labels if s not in low_quality]

    print(f"Passed:    {len(high_quality)}/{n}")
    print(f"Discarded: {len(low_quality)}/{n}")

    # Write filtered matrix
    hq_set = set(high_quality)
    hq_idx = [i for i, s in enumerate(labels) if s in hq_set]
    filtered = matrix[np.ix_(hq_idx, hq_idx)]

    with out_matrix.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([""] + high_quality)
        for i, name in enumerate(high_quality):
            writer.writerow([name] + [f"{v:.6f}" for v in filtered[i]])

    # Write low quality list
    with out_low_quality.open("w") as fh:
        fh.write("# SAGs removed: insufficient MinHash connections\n")
        fh.write(f"# min_similarity={cfg['min_similarity']}, "
                 f"min_connections={cfg['min_connections']}\n")
        for sag in low_quality:
            fh.write(f"{sag}\t{connectivity[sag]}\n")

    # Write report
    conn_vals = list(connectivity.values())
    with out_report.open("w") as fh:
        fh.write("# Sourmash Quality Filter Report\n")
        fh.write(f"min_similarity:            {cfg['min_similarity']}\n")
        fh.write(f"min_connections:           {cfg['min_connections']}\n")
        fh.write(f"total_sags:                {n}\n")
        fh.write(f"passed:                    {len(high_quality)}\n")
        fh.write(f"discarded:                 {len(low_quality)}\n")
        fh.write(f"retention_pct:             "
                 f"{len(high_quality)/n*100:.1f}\n")
        fh.write(f"mean_connections:          {np.mean(conn_vals):.2f}\n")
        fh.write(f"max_connections_observed:  {max(conn_vals)}\n")
        fh.write(f"min_connections_observed:  {min(conn_vals)}\n")

    if len(high_quality) == 0:
        print("ERROR: No SAGs passed quality filter. "
              "Consider lowering --min_similarity or --min_connections.")
        raise SystemExit(1)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--matrix_csv",      required=True, type=Path)
    p.add_argument("--out_matrix",      required=True, type=Path)
    p.add_argument("--out_report",      required=True, type=Path)
    p.add_argument("--out_low_quality", required=True, type=Path)
    p.add_argument("--min_similarity",  type=float,
                   default=DEFAULTS["min_similarity"])
    p.add_argument("--min_connections", type=int,
                   default=DEFAULTS["min_connections"])
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    cfg  = {k: getattr(args, k) for k in DEFAULTS}
    run(
        matrix_csv      = args.matrix_csv,
        out_matrix      = args.out_matrix,
        out_report      = args.out_report,
        out_low_quality = args.out_low_quality,
        cfg             = cfg,
    )