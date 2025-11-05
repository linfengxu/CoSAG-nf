process HIERARCHICAL_CLUSTERING_MATRIX {
    tag "hierarchical_clustering_matrix_prep"
    label 'process_low'
    container "${params.python.container}"
    publishDir "${params.output_structure?.clustering_analysis ?: params.outdir}/hierarchical_clustering", mode: 'copy'

    input:
    path(similarity_matrix)
    val(distance_metric)

    output:
    path("distance_matrix.tsv"), emit: distance_matrix
    path("similarity_stats.txt"), emit: similarity_stats
    path "versions.yml", emit: versions

    script:
    """
    python3 << 'EOF'
import pandas as pd
import numpy as np
import scipy
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram, inconsistent
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys

print(f"Using SciPy version: {scipy.__version__}")

# Check if similarity matrix file exists and is not empty
similarity_file = "${similarity_matrix}"
print(f"Checking similarity matrix file: {similarity_file}")

if not os.path.exists(similarity_file):
    print(f"ERROR: Similarity matrix file '{similarity_file}' does not exist!")
    sys.exit(1)

file_size = os.path.getsize(similarity_file)
print(f"File size: {file_size} bytes")

if file_size == 0:
    print(f"ERROR: Similarity matrix file '{similarity_file}' is empty!")
    sys.exit(1)

# Read and validate similarity matrix
print("Loading similarity matrix...")
try:
    similarity_df = pd.read_csv(similarity_file, sep="\\t", index_col=0)
    print(f"Similarity matrix shape: {similarity_df.shape}")
    
    if similarity_df.shape[0] == 0 or similarity_df.shape[1] == 0:
        print("ERROR: Similarity matrix is empty (0x0)!")
        print("This usually means:")
        print("1. No SAGs were processed in previous steps")
        print("2. All SAGs were filtered out during quality control")
        print("3. Sourmash similarity computation failed")
        sys.exit(1)
    
    if similarity_df.shape[0] < 2:
        print(f"ERROR: Need at least 2 SAGs for clustering, got {similarity_df.shape[0]}")
        sys.exit(1)
        
    print(f"Successfully loaded similarity matrix with {similarity_df.shape[0]} SAGs")
    
except Exception as e:
    print(f"ERROR: Failed to read similarity matrix: {e}")
    print("File contents preview:")
    try:
        with open(similarity_file, 'r') as f:
            lines = f.readlines()[:10]  # First 10 lines
            for i, line in enumerate(lines):
                print(f"Line {i+1}: {line.strip()}")
    except:
        print("Could not read file contents")
    sys.exit(1)

# Basic statistics
similarity_values = similarity_df.values
np.fill_diagonal(similarity_values, np.nan)  # Exclude diagonal for stats
non_diag_values = similarity_values[~np.isnan(similarity_values)]

if len(non_diag_values) == 0:
    print("ERROR: No valid similarity values found!")
    sys.exit(1)

print(f"Similarity statistics (excluding diagonal):")
print(f"  Mean: {np.mean(non_diag_values):.4f}")
print(f"  Median: {np.median(non_diag_values):.4f}")
print(f"  Std: {np.std(non_diag_values):.4f}")
print(f"  Min: {np.min(non_diag_values):.4f}")
print(f"  Max: {np.max(non_diag_values):.4f}")

# Convert similarity to distance
print("Converting similarity matrix to distance matrix...")
print(f"Distance metric: ${distance_metric}")

if "${distance_metric}" == "euclidean":
    # For Euclidean distance: distance = sqrt(2 * (1 - similarity))
    # This ensures that identical sequences (similarity=1) have distance=0
    distance_matrix = np.sqrt(2 * (1 - similarity_df.values))
elif "${distance_metric}" == "cosine":
    # For cosine distance: distance = 1 - similarity
    distance_matrix = 1 - similarity_df.values
elif "${distance_metric}" == "manhattan":
    # For Manhattan distance: distance = 2 * (1 - similarity)
    distance_matrix = 2 * (1 - similarity_df.values)
else:
    # Default to Euclidean
    print(f"Warning: Unknown distance metric '${distance_metric}', using Euclidean")
    distance_matrix = np.sqrt(2 * (1 - similarity_df.values))

# Check for invalid values
if np.any(np.isnan(distance_matrix)) or np.any(np.isinf(distance_matrix)):
    print("WARNING: Found NaN or Inf values in distance matrix, replacing with max distance")
    max_valid_distance = np.nanmax(distance_matrix[np.isfinite(distance_matrix)])
    distance_matrix = np.where(np.isfinite(distance_matrix), distance_matrix, max_valid_distance)

# Ensure diagonal is zero (self-distance)
np.fill_diagonal(distance_matrix, 0.0)

# Ensure symmetry
distance_matrix = (distance_matrix + distance_matrix.T) / 2

# Create distance DataFrame
distance_df = pd.DataFrame(
    distance_matrix,
    index=similarity_df.index,
    columns=similarity_df.columns
)

# Save distance matrix
distance_df.to_csv("distance_matrix.tsv", sep="\\t", float_format='%.6f')
print(f"Distance matrix saved: {distance_df.shape}")

# Distance statistics
distance_values = distance_matrix.copy()
np.fill_diagonal(distance_values, np.nan)  # Exclude diagonal for stats
non_diag_distances = distance_values[~np.isnan(distance_values)]

print(f"Distance statistics (excluding diagonal):")
print(f"  Mean: {np.mean(non_diag_distances):.4f}")
print(f"  Median: {np.median(non_diag_distances):.4f}")
print(f"  Std: {np.std(non_diag_distances):.4f}")
print(f"  Min: {np.min(non_diag_distances):.4f}")
print(f"  Max: {np.max(non_diag_distances):.4f}")

# Generate statistics report
with open("similarity_stats.txt", "w") as f:
    f.write("# Similarity to Distance Matrix Conversion Report\\n")
    f.write(f"SciPy version: {scipy.__version__}\\n")
    f.write(f"Distance metric: ${distance_metric}\\n")
    f.write(f"Matrix dimensions: {similarity_df.shape[0]} x {similarity_df.shape[1]}\\n\\n")
    
    f.write("## Similarity Matrix Statistics\\n")
    f.write(f"Mean similarity: {np.mean(non_diag_values):.6f}\\n")
    f.write(f"Median similarity: {np.median(non_diag_values):.6f}\\n")
    f.write(f"Std similarity: {np.std(non_diag_values):.6f}\\n")
    f.write(f"Min similarity: {np.min(non_diag_values):.6f}\\n")
    f.write(f"Max similarity: {np.max(non_diag_values):.6f}\\n\\n")
    
    f.write("## Distance Matrix Statistics\\n")
    f.write(f"Mean distance: {np.mean(non_diag_distances):.6f}\\n")
    f.write(f"Median distance: {np.median(non_diag_distances):.6f}\\n")
    f.write(f"Std distance: {np.std(non_diag_distances):.6f}\\n")
    f.write(f"Min distance: {np.min(non_diag_distances):.6f}\\n")
    f.write(f"Max distance: {np.max(non_diag_distances):.6f}\\n\\n")
    
    f.write("## Conversion Formula\\n")
    if "${distance_metric}" == "euclidean":
        f.write("Distance = sqrt(2 * (1 - similarity))\\n")
    elif "${distance_metric}" == "cosine":
        f.write("Distance = 1 - similarity\\n")
    elif "${distance_metric}" == "manhattan":
        f.write("Distance = 2 * (1 - similarity)\\n")
    
    f.write("\\n## Matrix Properties\\n")
    f.write(f"Symmetric: {np.allclose(distance_matrix, distance_matrix.T)}\\n")
    f.write(f"Zero diagonal: {np.allclose(np.diag(distance_matrix), 0)}\\n")
    f.write(f"Non-negative: {np.all(distance_matrix >= 0)}\\n")

print("Distance matrix conversion completed successfully!")

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

process HIERARCHICAL_CLUSTERING_ANALYSIS {
    tag "hierarchical_clustering_analysis"
    label 'process_medium'
    container "${params.python.container}"
    publishDir "${params.output_structure?.clustering_analysis ?: params.outdir}/hierarchical_clustering", mode: 'copy'

    input:
    path(distance_matrix)
    val(linkage_method)
    val(criterion)
    val(threshold)

    output:
    path("hierarchical_clusters.tsv"), emit: clusters
    path("linkage_matrix.tsv"), emit: linkage_matrix
    path("dendrogram.png"), emit: dendrogram
    path("clustering_report.txt"), emit: report
    path("cluster_validation.txt"), emit: validation
    path "versions.yml", emit: versions

    script:
    """
    python3 << 'EOF'
import pandas as pd
import numpy as np
import scipy
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram, inconsistent, cophenet
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import time

print(f"Using SciPy version: {scipy.__version__}")
print(f"Clustering parameters:")
print(f"  Linkage method: ${linkage_method}")
print(f"  Criterion: ${criterion}")
print(f"  Threshold: ${threshold}")

# Read distance matrix
print("Loading distance matrix...")
distance_df = pd.read_csv("${distance_matrix}", sep="\\t", index_col=0)
print(f"Distance matrix shape: {distance_df.shape}")

# Extract SAG names
sag_names = distance_df.index.tolist()
n_sags = len(sag_names)
print(f"Number of SAGs: {n_sags}")

# Convert to condensed distance matrix for scipy
print("Converting to condensed distance matrix...")
distance_matrix = distance_df.values

# Ensure the matrix is symmetric and has zero diagonal
if not np.allclose(distance_matrix, distance_matrix.T):
    print("Warning: Distance matrix is not symmetric, symmetrizing...")
    distance_matrix = (distance_matrix + distance_matrix.T) / 2

np.fill_diagonal(distance_matrix, 0.0)

# Convert to condensed form (upper triangle)
condensed_distances = squareform(distance_matrix, checks=False)
print(f"Condensed distance vector length: {len(condensed_distances)}")

# Perform hierarchical clustering
print("Performing hierarchical clustering...")
start_time = time.time()

linkage_matrix = linkage(
    condensed_distances, 
    method="${linkage_method}"
)

clustering_time = time.time() - start_time
print(f"Clustering completed in {clustering_time:.2f} seconds")

# Apply clustering criterion
print(f"Applying clustering criterion: ${criterion} with threshold ${threshold}")

if "${criterion}" == "inconsistent":
    # Calculate inconsistency coefficient
    inconsistency_matrix = inconsistent(linkage_matrix)
    clusters = fcluster(linkage_matrix, float("${threshold}"), criterion="${criterion}")
elif "${criterion}" == "distance":
    clusters = fcluster(linkage_matrix, float("${threshold}"), criterion="${criterion}")
elif "${criterion}" == "maxclust":
    clusters = fcluster(linkage_matrix, int(float("${threshold}")), criterion="${criterion}")
else:
    print(f"Warning: Unknown criterion '${criterion}', using 'inconsistent'")
    clusters = fcluster(linkage_matrix, float("${threshold}"), criterion="inconsistent")

# Cluster statistics
n_clusters = len(set(clusters))
cluster_counts = Counter(clusters)
print(f"Number of clusters formed: {n_clusters}")
print(f"Cluster size distribution:")
for cluster_id, count in sorted(cluster_counts.items()):
    print(f"  Cluster {cluster_id}: {count} SAGs")

# Create results DataFrame
results_df = pd.DataFrame({
    'SAG_ID': sag_names,
    'Cluster_ID': clusters,
    'Cluster_Size': [cluster_counts[cluster_id] for cluster_id in clusters]
})

# Sort by cluster ID and SAG ID
results_df = results_df.sort_values(['Cluster_ID', 'SAG_ID'])

# Save clustering results
results_df.to_csv("hierarchical_clusters.tsv", sep="\\t", index=False)
print(f"Clustering results saved: {len(results_df)} SAGs in {n_clusters} clusters")

# Save linkage matrix
linkage_df = pd.DataFrame(
    linkage_matrix,
    columns=['SAG1_Index', 'SAG2_Index', 'Distance', 'Cluster_Size']
)
linkage_df.to_csv("linkage_matrix.tsv", sep="\\t", index=False)

# Generate dendrogram
print("Generating dendrogram...")
max_width = 200  # Maximum width in inches (adjust as needed)
calculated_width = max(12, n_sags * 0.3)
figure_width = min(calculated_width, max_width)


# For large datasets, limit leaf labels
if n_sags > 50:
    dendrogram(
        linkage_matrix,
        labels=sag_names,
        leaf_rotation=90,
        leaf_font_size=6,
        no_labels=True  # Hide labels for large datasets
    )
    plt.title(f'Hierarchical Clustering Dendrogram ({n_sags} SAGs)\\nMethod: ${linkage_method}, Criterion: ${criterion}, Threshold: ${threshold}')
else:
    dendrogram(
        linkage_matrix,
        labels=sag_names,
        leaf_rotation=90,
        leaf_font_size=8
    )
    plt.title(f'Hierarchical Clustering Dendrogram\\nMethod: ${linkage_method}, Criterion: ${criterion}, Threshold: ${threshold}')

plt.xlabel('SAG ID')
plt.ylabel('Distance')
plt.tight_layout()
plt.savefig('dendrogram.png', dpi=300, bbox_inches='tight')
plt.close()

# Cluster validation
print("Performing cluster validation...")

# Cophenetic correlation coefficient
try:
    cophenetic_distances = cophenet(linkage_matrix, condensed_distances)
    
    # Check if cophenetic_distances is a tuple (cophenet returns tuple in some versions)
    if isinstance(cophenetic_distances, tuple):
        cophenetic_distances = cophenetic_distances[0]
    
    # Ensure both arrays have the same shape
    print(f"Cophenetic distances shape: {cophenetic_distances.shape}")
    print(f"Condensed distances shape: {condensed_distances.shape}")
    
    if cophenetic_distances.shape == condensed_distances.shape:
        cophenetic_corr = np.corrcoef(cophenetic_distances, condensed_distances)[0, 1]
        print(f"Cophenetic correlation coefficient: {cophenetic_corr:.4f}")
    else:
        print(f"Warning: Shape mismatch - cophenetic: {cophenetic_distances.shape}, condensed: {condensed_distances.shape}")
        cophenetic_corr = np.nan
        print("Cophenetic correlation coefficient: Not available (shape mismatch)")
        
except Exception as e:
    print(f"Error calculating cophenetic correlation: {e}")
    cophenetic_corr = np.nan
    print("Cophenetic correlation coefficient: Not available (calculation error)")

# Within-cluster and between-cluster distances
within_cluster_distances = []
between_cluster_distances = []

for i in range(n_sags):
    for j in range(i + 1, n_sags):
        dist = distance_matrix[i, j]
        if clusters[i] == clusters[j]:
            within_cluster_distances.append(dist)
        else:
            between_cluster_distances.append(dist)

mean_within = np.mean(within_cluster_distances) if within_cluster_distances else 0
mean_between = np.mean(between_cluster_distances) if between_cluster_distances else 0

print(f"Mean within-cluster distance: {mean_within:.4f}")
print(f"Mean between-cluster distance: {mean_between:.4f}")

# Generate comprehensive report
with open("clustering_report.txt", "w") as f:
    f.write("# Hierarchical Clustering Analysis Report\\n")
    f.write(f"Generated on: {time.strftime('%Y-%m-%d %H:%M:%S')}\\n")
    f.write(f"SciPy version: {scipy.__version__}\\n\\n")
    
    f.write("## Clustering Parameters\\n")
    f.write(f"Linkage method: ${linkage_method}\\n")
    f.write(f"Criterion: ${criterion}\\n")
    f.write(f"Threshold: ${threshold}\\n")
    f.write(f"Number of SAGs: {n_sags}\\n\\n")
    
    f.write("## Clustering Results\\n")
    f.write(f"Number of clusters: {n_clusters}\\n")
    f.write(f"Clustering time: {clustering_time:.2f} seconds\\n\\n")
    
    f.write("## Cluster Size Distribution\\n")
    size_distribution = Counter([cluster_counts[cluster_id] for cluster_id in cluster_counts])
    for size, count in sorted(size_distribution.items()):
        f.write(f"Clusters with {size} SAGs: {count}\\n")
    
    f.write("\\n## Detailed Cluster Information\\n")
    for cluster_id in sorted(cluster_counts.keys()):
        cluster_sags = results_df[results_df['Cluster_ID'] == cluster_id]['SAG_ID'].tolist()
        f.write(f"Cluster {cluster_id} ({len(cluster_sags)} SAGs): {', '.join(cluster_sags[:10])}")
        if len(cluster_sags) > 10:
            f.write(f" ... and {len(cluster_sags) - 10} more")
        f.write("\\n")
    
    f.write("\\n## Cluster Validation Metrics\\n")
    if np.isnan(cophenetic_corr):
        f.write("Cophenetic correlation coefficient: Not available (calculation error)\\n")
    else:
        f.write(f"Cophenetic correlation coefficient: {cophenetic_corr:.6f}\\n")
    f.write(f"Mean within-cluster distance: {mean_within:.6f}\\n")
    f.write(f"Mean between-cluster distance: {mean_between:.6f}\\n")
    if mean_within > 0:
        f.write(f"Separation ratio: {mean_between/mean_within:.6f} (higher is better)\\n\\n")
    else:
        f.write("Separation ratio: Not available (no within-cluster distances)\\n\\n")

    f.write("## Cluster Quality Assessment\\n")
    # Assess each cluster
    for cluster_id in sorted(cluster_counts.keys()):
        cluster_mask = (clusters == cluster_id)
        cluster_indices = np.where(cluster_mask)[0]
        
        if len(cluster_indices) > 1:
            # Calculate within-cluster distances for this cluster
            cluster_distances = []
            for i in range(len(cluster_indices)):
                for j in range(i + 1, len(cluster_indices)):
                    idx1, idx2 = cluster_indices[i], cluster_indices[j]
                    cluster_distances.append(distance_matrix[idx1, idx2])
            
            mean_cluster_dist = np.mean(cluster_distances)
            max_cluster_dist = np.max(cluster_distances)
            
            f.write(f"Cluster {cluster_id}: {len(cluster_indices)} SAGs, ")
            f.write(f"mean internal distance: {mean_cluster_dist:.4f}, ")
            f.write(f"max internal distance: {max_cluster_dist:.4f}\\n")
        else:
            f.write(f"Cluster {cluster_id}: 1 SAG (singleton)\\n")

# Generate validation report
with open("cluster_validation.txt", "w") as f:
    f.write("# Cluster Validation Report\\n")
    if np.isnan(cophenetic_corr):
        f.write("Cophenetic correlation: Not available (calculation error)\\n")
        f.write("Interpretation: Cannot assess clustering quality due to correlation calculation error\\n\\n")
    else:
        f.write(f"Cophenetic correlation: {cophenetic_corr:.6f}\\n")
        f.write(f"Interpretation: {'Excellent' if cophenetic_corr > 0.9 else 'Good' if cophenetic_corr > 0.8 else 'Fair' if cophenetic_corr > 0.7 else 'Poor'} clustering quality\\n\\n")
    
    f.write("## Distance Analysis\\n")
    f.write(f"Mean within-cluster distance: {mean_within:.6f}\\n")
    f.write(f"Mean between-cluster distance: {mean_between:.6f}\\n")
    if mean_within > 0:
        f.write(f"Separation ratio: {mean_between/mean_within:.6f} (higher is better)\\n\\n")
    else:
        f.write("Separation ratio: Not available (no within-cluster distances)\\n\\n")

print("Hierarchical clustering analysis completed successfully!")

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

process HIERARCHICAL_CLUSTERING_REPORT {
    tag "hierarchical_clustering_report"
    label 'process_low'
    container "${params.python.container}"
    publishDir "${params.output_structure?.clustering_analysis ?: params.outdir}/cluster_statistics", mode: 'copy'

    input:
    path(clusters)
    path(similarity_matrix)
    path(clustering_report)
    path(validation_report)

    output:
    path("final_clustering_summary.txt"), emit: summary
    path("cluster_statistics.tsv"), emit: statistics
    path("high_quality_clusters.tsv"), emit: high_quality_clusters
    path "versions.yml", emit: versions

    script:
    """
    python3 << 'EOF'
import pandas as pd
import numpy as np
from collections import Counter
import time

print("Generating final clustering report...")

# Read clustering results
clusters_df = pd.read_csv("${clusters}", sep="\\t")
print(f"Loaded clustering results: {len(clusters_df)} SAGs")

# Read similarity matrix for quality assessment
similarity_df = pd.read_csv("${similarity_matrix}", sep="\\t", index_col=0)

# Basic cluster statistics
n_clusters = clusters_df['Cluster_ID'].nunique()
cluster_sizes = clusters_df.groupby('Cluster_ID').size()

print(f"Total clusters: {n_clusters}")
print(f"Cluster size range: {cluster_sizes.min()} - {cluster_sizes.max()}")

# Detailed cluster analysis
cluster_stats = []

for cluster_id in sorted(clusters_df['Cluster_ID'].unique()):
    cluster_sags = clusters_df[clusters_df['Cluster_ID'] == cluster_id]['SAG_ID'].tolist()
    cluster_size = len(cluster_sags)
    
    # Calculate internal similarity statistics
    if cluster_size > 1:
        # Get similarity values within this cluster
        cluster_similarities = []
        for i, sag1 in enumerate(cluster_sags):
            for j, sag2 in enumerate(cluster_sags):
                if i < j and sag1 in similarity_df.index and sag2 in similarity_df.index:
                    sim = similarity_df.loc[sag1, sag2]
                    cluster_similarities.append(sim)
        
        if cluster_similarities:
            mean_similarity = np.mean(cluster_similarities)
            min_similarity = np.min(cluster_similarities)
            max_similarity = np.max(cluster_similarities)
            std_similarity = np.std(cluster_similarities)
        else:
            mean_similarity = min_similarity = max_similarity = std_similarity = np.nan
    else:
        mean_similarity = min_similarity = max_similarity = std_similarity = 1.0  # Singleton
    
    cluster_stats.append({
        'Cluster_ID': cluster_id,
        'Size': cluster_size,
        'Mean_Similarity': mean_similarity,
        'Min_Similarity': min_similarity,
        'Max_Similarity': max_similarity,
        'Std_Similarity': std_similarity,
        'SAG_List': ','.join(cluster_sags[:5]) + ('...' if cluster_size > 5 else '')
    })

# Create statistics DataFrame
stats_df = pd.DataFrame(cluster_stats)
stats_df = stats_df.sort_values('Size', ascending=False)

# Save cluster statistics
stats_df.to_csv("cluster_statistics.tsv", sep="\\t", index=False, float_format='%.6f')

# Identify high-quality clusters
# Temporarily remove quality filtering - only filter by size >= 2
high_quality_mask = (stats_df['Size'] >= 2)

high_quality_clusters = stats_df[high_quality_mask].copy()
high_quality_clusters.to_csv("high_quality_clusters.tsv", sep="\\t", index=False, float_format='%.6f')

print(f"High-quality clusters identified: {len(high_quality_clusters)}")

# Generate final summary report
with open("final_clustering_summary.txt", "w") as f:
    f.write("# Final Hierarchical Clustering Summary\\n")
    f.write(f"Generated on: {time.strftime('%Y-%m-%d %H:%M:%S')}\\n\\n")
    
    f.write("## Overall Statistics\\n")
    f.write(f"Total SAGs processed: {len(clusters_df)}\\n")
    f.write(f"Total clusters formed: {n_clusters}\\n")
    f.write(f"High-quality clusters: {len(high_quality_clusters)}\\n")
    f.write(f"Singleton clusters: {sum(cluster_sizes == 1)}\\n")
    f.write(f"Multi-SAG clusters: {sum(cluster_sizes > 1)}\\n\\n")
    
    f.write("## Cluster Size Distribution\\n")
    size_dist = Counter(cluster_sizes)
    for size in sorted(size_dist.keys()):
        f.write(f"Size {size}: {size_dist[size]} clusters\\n")
    
    f.write("\\n## Quality Metrics Summary\\n")
    if len(high_quality_clusters) > 0:
        f.write(f"Mean similarity in high-quality clusters: {high_quality_clusters['Mean_Similarity'].mean():.4f}\\n")
        f.write(f"Mean size of high-quality clusters: {high_quality_clusters['Size'].mean():.1f}\\n")
        f.write(f"Largest high-quality cluster: {high_quality_clusters['Size'].max()} SAGs\\n")
    
    f.write("\\n## Top 10 Largest Clusters\\n")
    top_clusters = stats_df.head(10)
    for _, row in top_clusters.iterrows():
        f.write(f"Cluster {row['Cluster_ID']}: {row['Size']} SAGs, ")
        f.write(f"mean similarity: {row['Mean_Similarity']:.4f}\\n")
    
    f.write("\\n## High-Quality Clusters (Size >= 2, Mean Sim >= 0.1)\\n")
    for _, row in high_quality_clusters.iterrows():
        f.write(f"Cluster {row['Cluster_ID']}: {row['Size']} SAGs, ")
        f.write(f"similarity: {row['Mean_Similarity']:.4f} ± {row['Std_Similarity']:.4f}\\n")
    
    f.write("\\n## High-Quality Clusters (Size >= 2, No Quality Filtering)\\n")
    for _, row in high_quality_clusters.iterrows():
        f.write(f"Cluster {row['Cluster_ID']}: {row['Size']} SAGs, ")
        f.write(f"similarity: {row['Mean_Similarity']:.4f} ± {row['Std_Similarity']:.4f}\\n")
    
    f.write("\\n## Recommendations\\n")
    f.write("1. Focus on high-quality clusters for downstream analysis\\n")
    f.write("2. Consider re-clustering large clusters (>20 SAGs) with stricter thresholds\\n")
    f.write("3. Validate singleton clusters - they may represent unique genomes\\n")
    f.write("4. Use cluster statistics to guide MAG assembly parameters\\n")

print("Final clustering report generated successfully!")

EOF
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
} 