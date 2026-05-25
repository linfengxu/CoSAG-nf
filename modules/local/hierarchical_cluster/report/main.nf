process HIERARCHICAL_CLUSTER_REPORT_NF {
    tag "cluster_report"
    label 'process_low'

   // container 'python:3.11-slim'   // 替换为你的 Python 镜像
    //container '/cpfs01/projects-HDD/cfff-86962b7a8e68_HDD/public/singularity_sif/python3.8_v1_Bio_updated.sif'
    container 'quay.io/xulf2022/python3.8_bio:v1'
    input:
    path(clusters_tsv)
    path(similarity_csv)

    output:
    path("final_clustering_summary.txt"),  emit: summary
    path("cluster_statistics.tsv"),        emit: statistics
    path("high_quality_clusters.tsv"),     emit: high_quality_clusters
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    hierarchical_cluster_report.py \\
        --clusters_tsv              ${clusters_tsv} \\
        --similarity_csv            ${similarity_csv} \\
        --out_summary               final_clustering_summary.txt \\
        --out_statistics            cluster_statistics.tsv \\
        --out_high_quality          high_quality_clusters.tsv \\
        --min_cluster_size          ${params.cluster_min_size} \\
        --max_cluster_size          ${params.cluster_max_size} \\
        --min_similarity_in_cluster ${params.cluster_min_similarity_in_cluster} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}

process HIERARCHICAL_CLUSTER_REPORT {
    tag "hierarchical_clustering_report"
    label 'process_low'
    container 'quay.io/xulf2022/python3.8_bio:v1'
    //publishDir "${params.output_structure?.clustering_analysis ?: params.outdir}/cluster_statistics", mode: 'copy'

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
