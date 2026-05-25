process HIERARCHICAL_CLUSTERING_ANALYSIS {
    tag "hierarchical_clustering_analysis"
    label 'process_medium'
    container 'quay.io/xulf2022/python3.8_bio:v1'

    input:
    path(distance_matrix)
    val(linkage_method)
    val(criterion)
    val(threshold)

    output:
    path("hierarchical_clusters.tsv"), emit: clusters
    path("linkage_matrix.tsv"),        emit: linkage_matrix
    path("dendrogram.png"),            emit: dendrogram
    path("clustering_report.txt"),     emit: report
    path("cluster_validation.txt"),    emit: validation
    path("versions.yml"),              emit: versions

    script:
    """
    python3 << 'EOF'
import time
from collections import Counter

import numpy as np
import pandas as pd
import scipy
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram, inconsistent, cophenet

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

print(f"Using SciPy version: {scipy.__version__}")
print("Clustering parameters:")
print(f"  Linkage method: ${linkage_method}")
print(f"  Criterion:      ${criterion}")
print(f"  Threshold:      ${threshold}")

# ---------------------------------------------------------------------------
# 1. Load distance matrix
# ---------------------------------------------------------------------------
distance_df = pd.read_csv("${distance_matrix}", sep="\\t", index_col=0)
print(f"Distance matrix shape: {distance_df.shape}")
sag_names = distance_df.index.tolist()
n_sags = len(sag_names)
print(f"Number of SAGs: {n_sags}")

D = distance_df.values.astype(np.float64)

# Enforce matrix properties
if not np.allclose(D, D.T):
    print("Warning: Distance matrix is not symmetric; symmetrizing.")
    D = (D + D.T) / 2.0
np.fill_diagonal(D, 0.0)

# ---------------------------------------------------------------------------
# 2. Sanity check: flag unexpected distance ranges
#    (catches formula regressions upstream, e.g. sqrt(2(1-J)) with max ~1.414)
# ---------------------------------------------------------------------------
off_diag = D[~np.eye(n_sags, dtype=bool)]
d_min, d_max = float(off_diag.min()), float(off_diag.max())
d_mean, d_std = float(off_diag.mean()), float(off_diag.std())
print(f"Distance range: [{d_min:.4f}, {d_max:.4f}] "
      f"(mean = {d_mean:.4f}, std = {d_std:.4f})")
if d_max > 1.0 + 1e-6:
    print(f"WARNING: Max distance {d_max:.4f} > 1.0. "
          "Expected <=1.0 for Jaccard-derived distances (1-J or 1-J^(1/k)). "
          "Check upstream distance formulation.")
if d_min < -1e-6:
    print(f"WARNING: Negative distance {d_min:.4f} detected.")

condensed = squareform(D, checks=False)
print(f"Condensed distance vector length: {len(condensed)}")

# ---------------------------------------------------------------------------
# 3. Hierarchical clustering
#    NOTE: scipy.cluster.hierarchy.linkage ignores the 'metric' argument when
#    given a 1D condensed distance vector; the distance values themselves
#    fully determine the clustering (together with the linkage method).
# ---------------------------------------------------------------------------
print("Performing hierarchical clustering...")
t0 = time.time()
linkage_matrix = linkage(condensed, method="${linkage_method}")
clustering_time = time.time() - t0
print(f"Clustering completed in {clustering_time:.2f} s")

# ---------------------------------------------------------------------------
# 4. Apply clustering criterion
# ---------------------------------------------------------------------------
criterion = "${criterion}".strip().lower()
threshold_str = "${threshold}"
print(f"Applying criterion '{criterion}' with threshold {threshold_str}")

if criterion == "inconsistent":
    _ = inconsistent(linkage_matrix)  # kept for side-effect / sanity
    clusters = fcluster(linkage_matrix, float(threshold_str), criterion="inconsistent")
elif criterion == "distance":
    clusters = fcluster(linkage_matrix, float(threshold_str), criterion="distance")
elif criterion == "maxclust":
    clusters = fcluster(linkage_matrix, int(float(threshold_str)), criterion="maxclust")
else:
    print(f"Warning: Unknown criterion '{criterion}', defaulting to 'inconsistent'.")
    clusters = fcluster(linkage_matrix, float(threshold_str), criterion="inconsistent")

# ---------------------------------------------------------------------------
# 5. Cluster statistics and outputs
# ---------------------------------------------------------------------------
cluster_counts = Counter(clusters)
n_clusters = len(cluster_counts)
print(f"Number of clusters: {n_clusters}")

results_df = pd.DataFrame({
    'SAG_ID':       sag_names,
    'Cluster_ID':   clusters,
    'Cluster_Size': [cluster_counts[c] for c in clusters],
}).sort_values(['Cluster_ID', 'SAG_ID'])
results_df.to_csv("hierarchical_clusters.tsv", sep="\\t", index=False)

pd.DataFrame(
    linkage_matrix,
    columns=['Node1_Index', 'Node2_Index', 'Distance', 'Cluster_Size'],
).to_csv("linkage_matrix.tsv", sep="\\t", index=False)

# ---------------------------------------------------------------------------
# 6. Dendrogram
# ---------------------------------------------------------------------------
print("Generating dendrogram...")
figure_width = min(max(12, n_sags * 0.3), 200)
plt.figure(figsize=(figure_width, 8))
if n_sags > 50:
    dendrogram(linkage_matrix, labels=sag_names,
               leaf_rotation=90, leaf_font_size=6, no_labels=True)
else:
    dendrogram(linkage_matrix, labels=sag_names,
               leaf_rotation=90, leaf_font_size=8)
plt.title(f'Hierarchical Clustering Dendrogram ({n_sags} SAGs)\\n'
          f'Method: ${linkage_method}, Criterion: {criterion}, Threshold: {threshold_str}')
plt.xlabel('SAG ID')
plt.ylabel('Distance')
plt.tight_layout()
plt.savefig('dendrogram.png', dpi=300, bbox_inches='tight')
plt.close()

# ---------------------------------------------------------------------------
# 7. Validation: cophenetic correlation + within/between-cluster distances
# ---------------------------------------------------------------------------
cophenetic_corr = np.nan
try:
    coph_out = cophenet(linkage_matrix, condensed)
    # cophenet can return either a float (corr) or a tuple (corr, coph_dists)
    if isinstance(coph_out, tuple):
        if len(coph_out) >= 2 and hasattr(coph_out[1], 'shape'):
            coph_dists = coph_out[1]
            if coph_dists.shape == condensed.shape:
                cophenetic_corr = float(np.corrcoef(coph_dists, condensed)[0, 1])
            else:
                print(f"Warning: shape mismatch coph={coph_dists.shape}, cond={condensed.shape}")
        else:
            cophenetic_corr = float(coph_out[0])
    else:
        cophenetic_corr = float(coph_out)
    print(f"Cophenetic correlation: {cophenetic_corr:.4f}")
except Exception as e:
    print(f"Error computing cophenetic correlation: {e}")

# Within / between cluster distances (vectorized for speed on large N)
cluster_arr = np.asarray(clusters)
iu = np.triu_indices(n_sags, k=1)
pair_d     = D[iu]
same_clust = (cluster_arr[iu[0]] == cluster_arr[iu[1]])

within  = pair_d[same_clust]
between = pair_d[~same_clust]
mean_within  = float(within.mean())  if within.size  else 0.0
mean_between = float(between.mean()) if between.size else 0.0
print(f"Mean within-cluster distance:  {mean_within:.4f}")
print(f"Mean between-cluster distance: {mean_between:.4f}")

# ---------------------------------------------------------------------------
# 8. Reports
# ---------------------------------------------------------------------------
with open("clustering_report.txt", "w") as f:
    f.write("# Hierarchical Clustering Analysis Report\\n")
    f.write(f"Generated on: {time.strftime('%Y-%m-%d %H:%M:%S')}\\n")
    f.write(f"SciPy version: {scipy.__version__}\\n\\n")

    f.write("## Clustering Parameters\\n")
    f.write(f"Linkage method: ${linkage_method}\\n")
    f.write(f"Criterion: {criterion}\\n")
    f.write(f"Threshold: {threshold_str}\\n")
    f.write(f"Number of SAGs: {n_sags}\\n\\n")

    f.write("## Input Distance Summary\\n")
    f.write(f"Range (off-diagonal): [{d_min:.6f}, {d_max:.6f}]\\n")
    f.write(f"Mean: {d_mean:.6f}\\n")
    f.write(f"Std:  {d_std:.6f}\\n\\n")

    f.write("## Clustering Results\\n")
    f.write(f"Number of clusters: {n_clusters}\\n")
    f.write(f"Clustering time: {clustering_time:.2f} s\\n\\n")

    f.write("## Cluster Size Distribution\\n")
    size_dist = Counter([cluster_counts[c] for c in cluster_counts])
    for size, count in sorted(size_dist.items()):
        f.write(f"Clusters with {size} SAGs: {count}\\n")

    f.write("\\n## Detailed Cluster Information\\n")
    for cluster_id in sorted(cluster_counts.keys()):
        members = results_df[results_df['Cluster_ID'] == cluster_id]['SAG_ID'].tolist()
        preview = ', '.join(members[:10])
        f.write(f"Cluster {cluster_id} ({len(members)} SAGs): {preview}")
        if len(members) > 10:
            f.write(f" ... and {len(members) - 10} more")
        f.write("\\n")

    f.write("\\n## Cluster Validation Metrics\\n")
    if np.isnan(cophenetic_corr):
        f.write("Cophenetic correlation coefficient: Not available\\n")
    else:
        f.write(f"Cophenetic correlation coefficient: {cophenetic_corr:.6f}\\n")
    f.write(f"Mean within-cluster distance:  {mean_within:.6f}\\n")
    f.write(f"Mean between-cluster distance: {mean_between:.6f}\\n")
    if mean_within > 0:
        f.write(f"Separation ratio (between/within): {mean_between/mean_within:.6f}\\n\\n")
    else:
        f.write("Separation ratio: Not available\\n\\n")

    f.write("## Per-Cluster Quality\\n")
    for cluster_id in sorted(cluster_counts.keys()):
        idxs = np.where(cluster_arr == cluster_id)[0]
        if len(idxs) > 1:
            sub = D[np.ix_(idxs, idxs)]
            iu2 = np.triu_indices(len(idxs), k=1)
            cvals = sub[iu2]
            f.write(f"Cluster {cluster_id}: {len(idxs)} SAGs, "
                    f"mean intra dist: {float(cvals.mean()):.4f}, "
                    f"max intra dist: {float(cvals.max()):.4f}\\n")
        else:
            f.write(f"Cluster {cluster_id}: 1 SAG (singleton)\\n")

with open("cluster_validation.txt", "w") as f:
    f.write("# Cluster Validation Report\\n")
    if np.isnan(cophenetic_corr):
        f.write("Cophenetic correlation: Not available\\n")
        f.write("Interpretation: cannot assess clustering quality\\n\\n")
    else:
        f.write(f"Cophenetic correlation: {cophenetic_corr:.6f}\\n")
        tier = ('Excellent' if cophenetic_corr > 0.9
                else 'Good'   if cophenetic_corr > 0.8
                else 'Fair'   if cophenetic_corr > 0.7
                else 'Poor')
        f.write(f"Interpretation: {tier} clustering quality\\n\\n")

    f.write("## Distance Analysis\\n")
    f.write(f"Mean within-cluster distance:  {mean_within:.6f}\\n")
    f.write(f"Mean between-cluster distance: {mean_between:.6f}\\n")
    if mean_within > 0:
        f.write(f"Separation ratio: {mean_between/mean_within:.6f} (higher is better)\\n")
    else:
        f.write("Separation ratio: Not available\\n")

print("Hierarchical clustering analysis completed successfully.")

EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        scipy: \$(python -c "import scipy; print(scipy.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """
}
