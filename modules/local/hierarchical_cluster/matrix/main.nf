// ============================================================================
// Module: Hierarchical Clustering (CoSAG-nf Module 3)
// Author: Linfeng Xu
// Purpose: Convert sourmash Jaccard similarity -> principled distance,
//          then perform hierarchical clustering.
//
// Distance formulation:
//   - Default: Mash-style ANI distance  D = 1 - J^(1/k)   [Ondov et al. 2016]
//   - Option : Standard Jaccard distance D = 1 - J
//
// Rationale: Jaccard similarity does not inhabit a Euclidean vector space,
// so transformations such as sqrt(2(1-J)) (chord/cosine geometry) are not
// theoretically justified. The Mash ANI distance provides an ANI-calibrated
// metric that is particularly robust for single-cell genomes with variable
// completeness and is scale-matched to downstream FastANI validation.
// ============================================================================


process HIERARCHICAL_CLUSTER_MATRIX {
    tag "cluster_matrix"
    label 'process_low'
    container 'quay.io/xulf2022/python3.8_bio:v1'

    input:
    path(similarity_matrix)
    val(distance_metric)   // 'ani' (default, recommended) | 'jaccard' | legacy: 'euclidean'/'cosine'/'manhattan'
    val(kmer_size)         // k-mer size used upstream in sourmash (must match)

    output:
    path("distance_matrix.tsv"),  emit: distance_matrix
    path("similarity_stats.txt"), emit: similarity_stats
    path("versions.yml"),         emit: versions

    script:
    """
    python3 << 'EOF'
import os
import sys
import numpy as np
import pandas as pd
import scipy

print(f"Using SciPy version: {scipy.__version__}")

# ---------------------------------------------------------------------------
# 1. Load and validate similarity matrix
# ---------------------------------------------------------------------------
similarity_file = "${similarity_matrix}"
print(f"Similarity matrix file: {similarity_file}")

if not os.path.exists(similarity_file):
    print(f"ERROR: '{similarity_file}' does not exist!")
    sys.exit(1)

file_size = os.path.getsize(similarity_file)
print(f"File size: {file_size} bytes")
if file_size == 0:
    print(f"ERROR: '{similarity_file}' is empty!")
    sys.exit(1)

try:
    similarity_df = pd.read_csv(similarity_file, sep="\\t", index_col=0)
except Exception as e:
    print(f"ERROR: Failed to read similarity matrix: {e}")
    try:
        with open(similarity_file) as f:
            for i, line in enumerate(f.readlines()[:10]):
                print(f"Line {i+1}: {line.strip()}")
    except Exception:
        pass
    sys.exit(1)

print(f"Similarity matrix shape: {similarity_df.shape}")
if similarity_df.shape[0] < 2 or similarity_df.shape[0] != similarity_df.shape[1]:
    print(f"ERROR: Need a square matrix with >=2 SAGs, got {similarity_df.shape}")
    sys.exit(1)

S = similarity_df.values.astype(np.float64)

# ---------------------------------------------------------------------------
# 2. Similarity statistics
# ---------------------------------------------------------------------------
S_off = S.copy()
np.fill_diagonal(S_off, np.nan)
s_vals = S_off[~np.isnan(S_off)]
if len(s_vals) == 0:
    print("ERROR: No valid off-diagonal similarity values!")
    sys.exit(1)

print("Similarity statistics (off-diagonal):")
print(f"  mean   = {np.mean(s_vals):.6f}")
print(f"  median = {np.median(s_vals):.6f}")
print(f"  std    = {np.std(s_vals):.6f}")
print(f"  min    = {np.min(s_vals):.6f}")
print(f"  max    = {np.max(s_vals):.6f}")

# Clamp similarities to [0,1] to guard against tiny numerical excursions
S_clipped = np.clip(S, 0.0, 1.0)

# ---------------------------------------------------------------------------
# 3. Convert similarity -> distance
# ---------------------------------------------------------------------------
requested_metric = "${distance_metric}".strip().lower()
kmer_size = int("${kmer_size}")
if kmer_size < 1:
    print(f"ERROR: Invalid kmer_size={kmer_size}; must be >= 1.")
    sys.exit(1)

# Map legacy metric names with a deprecation warning
legacy_map = {
    "euclidean":  ("ani",     "Legacy 'euclidean' used sqrt(2(1-J)), which lacks theoretical basis for Jaccard similarity. Mapping to 'ani' (recommended)."),
    "cosine":     ("jaccard", "Legacy 'cosine' computed 1-J, which is actually the Jaccard distance. Mapping to 'jaccard'."),
    "manhattan":  ("jaccard", "Legacy 'manhattan' computed 2(1-J), a rescaled Jaccard distance. Mapping to 'jaccard'."),
}
effective_metric = requested_metric
if requested_metric in legacy_map:
    new_metric, msg = legacy_map[requested_metric]
    print(f"WARNING (deprecation): distance_metric='{requested_metric}' is deprecated.")
    print(f"WARNING: {msg}")
    effective_metric = new_metric

if effective_metric == "ani":
    # Mash ANI distance: D = 1 - J^(1/k)  (Ondov et al. 2016, Genome Biology)
    with np.errstate(invalid='ignore'):
        ani = np.power(S_clipped, 1.0 / kmer_size)
    D = 1.0 - ani
    formula_str = f"D = 1 - J^(1/k)  with k={kmer_size}  [Mash ANI distance, Ondov et al. 2016]"
elif effective_metric == "jaccard":
    # Standard Jaccard distance: D = 1 - J  (a true metric)
    D = 1.0 - S_clipped
    formula_str = "D = 1 - J  [Standard Jaccard distance]"
else:
    print(f"ERROR: Unsupported distance_metric '{requested_metric}'.")
    print("       Supported: 'ani' (recommended), 'jaccard'.")
    print("       Deprecated (auto-mapped): 'euclidean','cosine','manhattan'.")
    sys.exit(1)

print(f"Distance formula: {formula_str}")

# ---------------------------------------------------------------------------
# 4. Sanity checks and cleanup
# ---------------------------------------------------------------------------
if not np.all(np.isfinite(D)):
    print("WARNING: Non-finite values in distance matrix; replacing with max finite distance.")
    max_finite = np.nanmax(D[np.isfinite(D)]) if np.any(np.isfinite(D)) else 1.0
    D = np.where(np.isfinite(D), D, max_finite)

D = np.clip(D, 0.0, None)            # enforce non-negativity
np.fill_diagonal(D, 0.0)              # self-distance = 0
D = (D + D.T) / 2.0                   # enforce symmetry

# Range sanity check (helps catch formula regressions in the future)
d_max = float(D[~np.eye(D.shape[0], dtype=bool)].max())
if d_max > 1.0 + 1e-6:
    print(f"WARNING: Max off-diagonal distance = {d_max:.4f} > 1.0, "
          "which is unusual for Jaccard-derived distances. "
          "Check upstream similarity matrix and metric choice.")

# ---------------------------------------------------------------------------
# 5. Save outputs
# ---------------------------------------------------------------------------
distance_df = pd.DataFrame(D, index=similarity_df.index, columns=similarity_df.columns)
distance_df.to_csv("distance_matrix.tsv", sep="\\t", float_format="%.6f")
print(f"Distance matrix saved: {distance_df.shape}")

D_off = D.copy()
np.fill_diagonal(D_off, np.nan)
d_vals = D_off[~np.isnan(D_off)]

print("Distance statistics (off-diagonal):")
print(f"  mean   = {np.mean(d_vals):.6f}")
print(f"  median = {np.median(d_vals):.6f}")
print(f"  std    = {np.std(d_vals):.6f}")
print(f"  min    = {np.min(d_vals):.6f}")
print(f"  max    = {np.max(d_vals):.6f}")

with open("similarity_stats.txt", "w") as f:
    f.write("# Similarity-to-Distance Conversion Report\\n")
    f.write(f"SciPy version: {scipy.__version__}\\n")
    f.write(f"Requested metric: {requested_metric}\\n")
    f.write(f"Effective metric: {effective_metric}\\n")
    f.write(f"K-mer size (k): {kmer_size}\\n")
    f.write(f"Matrix dimensions: {similarity_df.shape[0]} x {similarity_df.shape[1]}\\n\\n")

    f.write("## Similarity Matrix Statistics (off-diagonal)\\n")
    f.write(f"Mean similarity:   {np.mean(s_vals):.6f}\\n")
    f.write(f"Median similarity: {np.median(s_vals):.6f}\\n")
    f.write(f"Std similarity:    {np.std(s_vals):.6f}\\n")
    f.write(f"Min similarity:    {np.min(s_vals):.6f}\\n")
    f.write(f"Max similarity:    {np.max(s_vals):.6f}\\n\\n")

    f.write("## Distance Matrix Statistics (off-diagonal)\\n")
    f.write(f"Mean distance:   {np.mean(d_vals):.6f}\\n")
    f.write(f"Median distance: {np.median(d_vals):.6f}\\n")
    f.write(f"Std distance:    {np.std(d_vals):.6f}\\n")
    f.write(f"Min distance:    {np.min(d_vals):.6f}\\n")
    f.write(f"Max distance:    {np.max(d_vals):.6f}\\n\\n")

    f.write("## Conversion Formula\\n")
    f.write(formula_str + "\\n\\n")

    f.write("## Matrix Properties\\n")
    f.write(f"Symmetric:    {np.allclose(D, D.T)}\\n")
    f.write(f"Zero diagonal:{np.allclose(np.diag(D), 0)}\\n")
    f.write(f"Non-negative: {bool(np.all(D >= 0))}\\n")
    f.write(f"Max off-diagonal distance: {d_max:.6f}\\n")

print("Similarity-to-distance conversion completed successfully.")

EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        scipy: \$(python -c "import scipy; print(scipy.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
