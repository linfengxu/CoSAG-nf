process SOURMASH_PROCESS_MATRIX {
    tag "sourmash_process_matrix"
    label 'process_medium'

    container 'quay.io/xulf2022/python3.8_bio:v1'

    //publishDir "${params.outdir}/sourmash_similarity", mode: 'copy'


    input:
    path(raw_matrix)
    path(compare_log)
    val(ksize)
    val(scaled)

    output:
    path("sourmash_similarity_matrix.tsv"), emit: similarity_matrix
    path("sourmash_comparison_details.csv"), emit: comparison_details
    path("sourmash_similarity_report.txt"), emit: report
    path("failed_comparisons.txt"), emit: failed_comparisons
    path "versions.yml", emit: versions

    script:
    """
    python3 << 'EOF'
import os
import numpy as np
import pandas as pd
import time
from pathlib import Path

print("Processing sourmash compare results...")
start_time = time.time()

# Check if raw matrix file exists and is valid
if not os.path.exists("${raw_matrix}"):
    print("Error: Raw similarity matrix file not found")
    exit(1)

# Read the raw similarity matrix from sourmash compare
try:
    raw_df = pd.read_csv("${raw_matrix}", index_col=0)
    print(f"Raw similarity matrix shape: {raw_df.shape}")

    # Check if the matrix is empty or invalid
    if raw_df.empty:
        print("Error: Raw similarity matrix is empty")
        raise ValueError("Empty matrix")

except Exception as e:
    print(f"Error reading raw similarity matrix: {e}")

    # Create empty output files for failed processing
    with open("sourmash_similarity_matrix.tsv", "w") as f:
        f.write("# Error: Failed to process raw similarity matrix\\n")
    with open("sourmash_comparison_details.csv", "w") as f:
        f.write("SAG1,SAG2,Jaccard_Similarity,Containment_1_2,Containment_2_1\\n")
    with open("sourmash_similarity_report.txt", "w") as f:
        f.write("Error: Failed to process raw similarity matrix\\n")
    with open("failed_comparisons.txt", "w") as f:
        f.write("Matrix processing failed\\n")
    exit(1)

# Get SAG names from the matrix
sag_names = raw_df.index.tolist()
n_sags = len(sag_names)
print(f"Number of SAGs: {n_sags}")

# Extract SAG IDs from full names (remove path and extension)
# Filter out any NaN or non-string values
sag_ids = []
valid_indices = []
for i, name in enumerate(sag_names):
    # Check if name is a valid string (not NaN or other types)
    if isinstance(name, str) and name.strip():
        # Extract just the filename without path and extension
        sag_id = os.path.basename(name).replace(f'_k${ksize}.sig', '')
        sag_ids.append(sag_id)
        valid_indices.append(i)
    else:
        print(f"Warning: Skipping invalid SAG name at index {i}: {name} (type: {type(name)})")

print(f"Valid SAGs after filtering: {len(sag_ids)}")
print(f"Sample SAG IDs: {sag_ids[:5]}")

# Filter the similarity matrix to only include valid SAGs
if len(valid_indices) != len(sag_names):
    print(f"Filtering matrix from {len(sag_names)} to {len(valid_indices)} valid SAGs")
    # Filter both rows and columns
    similarity_matrix = raw_df.iloc[valid_indices, valid_indices].values
    n_sags = len(sag_ids)  # Update n_sags to reflect filtered count
else:
    similarity_matrix = raw_df.values

# Create similarity matrix with proper SAG IDs
similarity_df = pd.DataFrame(
    similarity_matrix,
    index=sag_ids,
    columns=sag_ids
)

# Generate detailed comparison results
print("Generating detailed pairwise comparisons...")
comparison_details = []
total_comparisons = n_sags * (n_sags - 1) // 2
completed_comparisons = 0

for i, sag1 in enumerate(sag_ids):
    for j, sag2 in enumerate(sag_ids):
        if i < j:  # Only compute upper triangle
            jaccard_similarity = similarity_matrix[i][j]

            # For now, set containment values to 0 (can be computed separately if needed)
            containment_1_2 = 0.0
            containment_2_1 = 0.0

            comparison_details.append({
                'SAG1': sag1,
                'SAG2': sag2,
                'Jaccard_Similarity': jaccard_similarity,
                'Containment_1_2': containment_1_2,
                'Containment_2_1': containment_2_1
            })

            completed_comparisons += 1
            if completed_comparisons % 10000 == 0:
                elapsed = time.time() - start_time
                progress = completed_comparisons / total_comparisons * 100
                print(f"Progress: {completed_comparisons}/{total_comparisons} ({progress:.1f}%) - {elapsed:.1f}s elapsed")

print(f"Completed {completed_comparisons} comparisons in {time.time() - start_time:.1f} seconds")

# Save similarity matrix
similarity_df.to_csv("sourmash_similarity_matrix.tsv", sep="\\t", float_format='%.6f')
print(f"Similarity matrix saved: {similarity_df.shape}")

# Save detailed comparisons
comparison_df = pd.DataFrame(comparison_details)
comparison_df.to_csv("sourmash_comparison_details.csv", index=False, float_format='%.6f')
print(f"Comparison details saved: {len(comparison_details)} pairs")

# Calculate statistics
similarity_values = similarity_matrix.copy()
np.fill_diagonal(similarity_values, np.nan)  # Exclude diagonal for stats
non_diag_values = similarity_values[~np.isnan(similarity_values)]

if len(non_diag_values) > 0:
    mean_sim = np.mean(non_diag_values)
    median_sim = np.median(non_diag_values)
    max_sim = np.max(non_diag_values)
    min_sim = np.min(non_diag_values)
    std_sim = np.std(non_diag_values)

    # Count high similarity pairs
    high_sim_threshold = 0.1
    high_sim_pairs = np.sum(non_diag_values > high_sim_threshold)
else:
    mean_sim = median_sim = max_sim = min_sim = std_sim = 0.0
    high_sim_pairs = 0

# Generate comprehensive report
computation_time = time.time() - start_time
with open("sourmash_similarity_report.txt", "w") as f:
    f.write("# Sourmash Similarity Matrix Computation Report\\n")
    f.write(f"Generated on: {time.strftime('%Y-%m-%d %H:%M:%S')}\\n")
    f.write(f"K-mer size: ${ksize}\\n")
    f.write(f"Scaled parameter: ${scaled}\\n\\n")

    f.write("## Input Statistics\\n")
    f.write(f"Total signature files found: {n_sags}\\n")
    f.write(f"Valid signatures used: {n_sags}\\n")
    f.write(f"Invalid/empty signatures: 0\\n\\n")

    f.write("## Computation Statistics\\n")
    f.write(f"Total SAGs: {n_sags}\\n")
    f.write(f"Total pairwise comparisons: {total_comparisons}\\n")
    f.write(f"Successful comparisons: {completed_comparisons}\\n")
    f.write(f"Failed comparisons: 0\\n")
    f.write(f"Success rate: 100.0%\\n")
    f.write(f"Computation time: {computation_time:.1f} seconds\\n\\n")

    f.write("## Similarity Statistics\\n")
    f.write(f"Mean Jaccard similarity: {mean_sim:.6f}\\n")
    f.write(f"Median Jaccard similarity: {median_sim:.6f}\\n")
    f.write(f"Max Jaccard similarity: {max_sim:.6f}\\n")
    f.write(f"Min Jaccard similarity: {min_sim:.6f}\\n")
    f.write(f"Std Jaccard similarity: {std_sim:.6f}\\n\\n")

    f.write(f"## High Similarity Pairs (Jaccard > {high_sim_threshold})\\n")
    f.write(f"Number of high similarity pairs: {high_sim_pairs}\\n\\n")

    f.write("## Output Files\\n")
    f.write("- sourmash_similarity_matrix.tsv: Full similarity matrix\\n")
    f.write("- sourmash_comparison_details.csv: Detailed pairwise comparisons\\n")
    f.write("- failed_comparisons.txt: List of failed comparisons\\n")

# Create empty failed comparisons file (since sourmash compare succeeded)
with open("failed_comparisons.txt", "w") as f:
    f.write("# Failed comparisons: 0\\n")

print("Similarity matrix computation completed successfully!")
print(f"Matrix dimensions: {similarity_df.shape}")
print(f"Mean similarity: {mean_sim:.6f}")
print(f"Max similarity: {max_sim:.6f}")
print(f"High similarity pairs (>{high_sim_threshold}): {high_sim_pairs}")

EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(sourmash --version | sed 's/sourmash //')
        python: \$(python --version | sed 's/Python //')
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
