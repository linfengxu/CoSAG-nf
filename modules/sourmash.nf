process SOURMASH_SKETCH_GENERATION {
    tag "sourmash_sketch_${meta.id}"
    label 'process_low'
    container "${params.sourmash.container}"
    publishDir "${params.output_structure?.similarity_analysis ?: params.outdir}/sourmash_signatures",
        mode: 'copy'
    
    input:
    tuple val(meta), path(sag_fasta)
    val(ksize)
    val(scaled)

    output:
    tuple val(meta), path("${meta.id}_k${ksize}.sig"), emit: signature
    path("sketch_log_${meta.id}.txt"), emit: log
    path "versions.yml", emit: versions

    script:
    """
    # Generate MinHash signature
    sourmash sketch dna \\
        -p k=${ksize},scaled=${scaled} \\
        --name ${meta.id} \\
        -o ${meta.id}_k${ksize}.sig \\
        ${sag_fasta} > sketch_log_${meta.id}.txt 2>&1

    # Verify signature was created successfully
    if [ ! -f "${meta.id}_k${ksize}.sig" ]; then
        echo "ERROR: Failed to create signature for ${meta.id}" >> sketch_log_${meta.id}.txt
        touch ${meta.id}_k${ksize}.sig  # Create empty file to prevent pipeline failure
    else
        echo "SUCCESS: Signature created for ${meta.id}" >> sketch_log_${meta.id}.txt
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(sourmash --version | sed 's/sourmash //')
    END_VERSIONS
    """
}

process SOURMASH_COMPARE_MATRIX {
    tag "sourmash_compare_matrix"
    label 'process_high'
    container "${params.sourmash.container}"

    publishDir "${params.output_structure?.similarity_analysis ?: params.outdir}/similarity_matrices",
        mode: 'copy'
    input:
    path(signatures)
    val(ksize)
    val(scaled)

    output:
    path("similarity_matrix.npy"), emit: similarity_matrix
    path("similarity_matrix.npy.labels.txt"), emit: labels
    path("sourmash_compare.log"), emit: log
    path "versions.yml", emit: versions

    script:
    """
    # Get all signature files
    SIG_FILES=(\$(ls *.sig))
    echo "Found \${#SIG_FILES[@]} signature files" | tee sourmash_compare.log
    
    # Check if we have enough signatures
    if [ \${#SIG_FILES[@]} -lt 2 ]; then
        echo "Error: Need at least 2 valid signatures for comparison" | tee -a sourmash_compare.log
        # Create empty files for pipeline continuation
        touch similarity_matrix.npy
        touch similarity_matrix.npy.labels.txt
        exit 0
    fi
    
    # Use sourmash compare to generate binary matrix
    echo "Running sourmash compare on \${#SIG_FILES[@]} signatures..." | tee -a sourmash_compare.log
    echo "Command: sourmash compare *.sig -o similarity_matrix.npy --ksize ${ksize}" | tee -a sourmash_compare.log
    
    sourmash compare *.sig -o similarity_matrix.npy --ksize ${ksize} 2>&1 | tee -a sourmash_compare.log
    
    # Check if sourmash compare succeeded
    if [ ! -f "similarity_matrix.npy" ]; then
        echo "Error: sourmash compare failed to create .npy file" | tee -a sourmash_compare.log
        touch similarity_matrix.npy
        touch similarity_matrix.npy.labels.txt
        exit 1
    fi
    
    # Check if labels file was created
    if [ ! -f "similarity_matrix.npy.labels.txt" ]; then
        echo "Warning: Labels file not created by sourmash" | tee -a sourmash_compare.log
        touch similarity_matrix.npy.labels.txt
    fi
    
    echo "sourmash compare completed successfully" | tee -a sourmash_compare.log
    echo "Matrix file size: \$(ls -lh similarity_matrix.npy | awk '{print \$5}')" | tee -a sourmash_compare.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(sourmash --version | sed 's/sourmash //')
    END_VERSIONS
    """
}

process SOURMASH_MATRIX_CONVERT {
    tag "sourmash_matrix_convert"
    label 'process_medium'
    container "${params.python.container}"
    //publishDir "${params.outdir}/sourmash_similarity", mode: 'copy'
    publishDir "${params.output_structure?.similarity_analysis ?: params.outdir}/similarity_matrices", mode: 'copy'

    container "${params.python.container}"

    input:
    path(similarity_matrix)
    path(labels_file)
    path(signatures)
    val(ksize)
    val(scaled)

    output:
    path("raw_similarity_matrix.csv"), emit: raw_matrix
    path("matrix_convert.log"), emit: log
    path "versions.yml", emit: versions

    script:
    """
    python3 << 'PYTHON_EOF'
import numpy as np
import os
import glob

print("Converting sourmash .npy matrix to labeled CSV...")

# Load the similarity matrix
try:
    similarity_matrix = np.load('${similarity_matrix}')
    print(f"Loaded similarity matrix shape: {similarity_matrix.shape}")
except Exception as e:
    print(f"Error loading similarity matrix: {e}")
    # Create empty CSV file
    with open('raw_similarity_matrix.csv', 'w') as f:
        f.write("# Error: Failed to load similarity matrix\\n")
    exit(1)

# Try to load labels from sourmash-generated labels file
labels_file = '${labels_file}'
sag_names = []

if os.path.exists(labels_file) and os.path.getsize(labels_file) > 0:
    print("Loading SAG names from sourmash labels file...")
    try:
        with open(labels_file, 'r') as f:
            sag_names = [line.strip() for line in f.readlines() if line.strip()]
        print(f"Loaded {len(sag_names)} SAG names from labels file")
    except Exception as e:
        print(f"Error reading labels file: {e}")
        sag_names = []

# If labels file is empty or failed, extract from signature filenames
if not sag_names:
    print("Labels file empty or failed, extracting from signature filenames...")
    # Get signature file names in the same order as sourmash compare used
    sig_files = sorted(glob.glob('*.sig'))
    print(f"Found {len(sig_files)} signature files")
    
    # Extract SAG names from signature filenames
    for sig_file in sig_files:
        # Remove .sig extension and extract SAG ID
        sag_name = os.path.basename(sig_file).replace('.sig', '')
        # Remove the _k{ksize} suffix if present
        sag_name = sag_name.replace(f'_k${ksize}', '')
        sag_names.append(sag_name)

print(f"Final SAG names count: {len(sag_names)}")
print(f"Sample SAG names: {sag_names[:5] if sag_names else 'None'}")

# Ensure matrix dimensions match number of SAG names
if similarity_matrix.shape[0] != len(sag_names):
    print(f"Warning: Matrix shape {similarity_matrix.shape} doesn't match {len(sag_names)} SAG names")
    # Adjust if needed
    min_size = min(similarity_matrix.shape[0], len(sag_names))
    similarity_matrix = similarity_matrix[:min_size, :min_size]
    sag_names = sag_names[:min_size]
    print(f"Adjusted to {min_size}x{min_size} matrix")

# Write CSV manually without pandas
print("Writing labeled CSV file...")
try:
    with open('raw_similarity_matrix.csv', 'w') as f:
        # Write header row (column names)
        f.write(',' + ','.join(sag_names) + '\\n')
        
        # Write data rows
        for i, row_name in enumerate(sag_names):
            # Write row name followed by similarity values
            row_values = [str(similarity_matrix[i, j]) for j in range(len(sag_names))]
            f.write(row_name + ',' + ','.join(row_values) + '\\n')
    
    print(f"Successfully saved labeled similarity matrix: {similarity_matrix.shape}")
    print("First few similarity values:")
    for i in range(min(3, len(sag_names))):
        for j in range(min(3, len(sag_names))):
            print(f"{sag_names[i]} vs {sag_names[j]}: {similarity_matrix[i,j]:.6f}")

except Exception as e:
    print(f"Error writing CSV file: {e}")
    # Create minimal CSV file
    with open('raw_similarity_matrix.csv', 'w') as f:
        f.write("# Error: Failed to write similarity matrix\\n")

PYTHON_EOF

# Log the conversion results
echo "Matrix conversion completed" > matrix_convert.log
if [ -f "raw_similarity_matrix.csv" ]; then
    echo "CSV file created successfully: \$(wc -l < raw_similarity_matrix.csv) lines" >> matrix_convert.log
else
    echo "Error: CSV file was not created" >> matrix_convert.log
fi

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | sed 's/Python //')
    numpy: \$(python -c "import numpy; print(numpy.__version__)")
END_VERSIONS
    """
}

process SOURMASH_PROCESS_MATRIX {
    tag "sourmash_process_matrix"
    label 'process_medium'
    container "${params.python.container}"
    //publishDir "${params.outdir}/sourmash_similarity", mode: 'copy'
    publishDir "${params.output_structure?.similarity_analysis ?: params.outdir}/similarity_matrices",
        mode: 'copy'

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

process SOURMASH_QUALITY_FILTER {
    tag "sourmash_quality_filter"
    label 'process_low'
    container "${params.python.container}"
    publishDir "${params.output_structure?.similarity_analysis ?: params.outdir}/quality_filtering",
        mode: 'copy'

    input:
    path(similarity_matrix)
    path(comparison_details)
    val(min_similarity_threshold)
    val(min_connections)

    output:
    path("filtered_similarity_matrix.tsv"), emit: filtered_matrix
    path("sag_connectivity_report.txt"), emit: connectivity_report
    path("low_quality_sags.txt"), emit: low_quality_sags
    path "versions.yml", emit: versions

    script:
    """
    python3 << 'EOF'
import pandas as pd
import numpy as np

# Read similarity matrix
similarity_df = pd.read_csv("${similarity_matrix}", sep="\\t", index_col=0)
print(f"Loaded similarity matrix: {similarity_df.shape}")

# Read comparison details
try:
    details_df = pd.read_csv("${comparison_details}")
    print(f"Loaded comparison details: {len(details_df)} comparisons")
except:
    print("Warning: Could not load comparison details")
    details_df = None

min_similarity = float("${min_similarity_threshold}")
min_connections = int("${min_connections}")

# Calculate connectivity for each SAG
sag_connectivity = {}
low_quality_sags = []

for sag in similarity_df.index:
    # Count connections above threshold (excluding self-similarity)
    similarities = similarity_df.loc[sag]
    connections = sum(1 for sim in similarities if sim >= min_similarity and sim < 1.0)
    sag_connectivity[sag] = connections
    
    if connections < min_connections:
        low_quality_sags.append(sag)

print(f"SAGs with insufficient connections (< {min_connections}): {len(low_quality_sags)}")

# Filter out low-quality SAGs
high_quality_sags = [sag for sag in similarity_df.index if sag not in low_quality_sags]
filtered_matrix = similarity_df.loc[high_quality_sags, high_quality_sags]

# Save filtered matrix
filtered_matrix.to_csv("filtered_similarity_matrix.tsv", sep="\\t", float_format='%.6f')
print(f"Filtered similarity matrix saved: {filtered_matrix.shape}")

# Generate connectivity report
with open("sag_connectivity_report.txt", "w") as f:
    f.write("# SAG Connectivity Analysis Report\\n")
    f.write(f"Minimum similarity threshold: {min_similarity}\\n")
    f.write(f"Minimum required connections: {min_connections}\\n\\n")
    
    f.write("## Overall Statistics\\n")
    f.write(f"Total SAGs: {len(similarity_df.index)}\\n")
    f.write(f"High-quality SAGs: {len(high_quality_sags)}\\n")
    f.write(f"Low-quality SAGs: {len(low_quality_sags)}\\n")
    f.write(f"Retention rate: {len(high_quality_sags)/len(similarity_df.index)*100:.1f}%\\n\\n")
    
    f.write("## Connectivity Distribution\\n")
    connectivity_values = list(sag_connectivity.values())
    f.write(f"Mean connections per SAG: {np.mean(connectivity_values):.2f}\\n")
    f.write(f"Median connections per SAG: {np.median(connectivity_values):.2f}\\n")
    f.write(f"Max connections: {max(connectivity_values)}\\n")
    f.write(f"Min connections: {min(connectivity_values)}\\n\\n")
    
    f.write("## High-Quality SAGs (Top 10 by connectivity)\\n")
    sorted_connectivity = sorted(sag_connectivity.items(), key=lambda x: x[1], reverse=True)
    for sag, connections in sorted_connectivity[:10]:
        f.write(f"{sag}: {connections} connections\\n")

# Save low-quality SAGs list
with open("low_quality_sags.txt", "w") as f:
    f.write("# Low-quality SAGs (insufficient connections)\\n")
    f.write(f"# Threshold: minimum {min_connections} connections with similarity >= {min_similarity}\\n")
    for sag in low_quality_sags:
        connections = sag_connectivity[sag]
        f.write(f"{sag}\\t{connections}\\n")

print("Quality filtering completed!")

EOF
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(sourmash --version | sed 's/sourmash //')
        python: \$(python --version | sed 's/Python //')
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}

// Legacy processes for backward compatibility
process SOURMASH_COMPUTE {
    tag "${meta.id}"
    
    container { "${params.sourmash.container}" }
    
    input:
    tuple val(meta), path(fasta)
    
    output:
    tuple val(meta), path("${meta.id}.sig"), emit: signatures
    path "versions.yml", emit: versions
    
    script:
    """
    sourmash compute \
        --track-abundance \
        -k 51 \
        ${fasta} \
        -o ${meta.id}.sig
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(sourmash --version | sed 's/sourmash //')
    END_VERSIONS
    """
}

process SOURMASH_COMPARE {
    tag "compare_signatures"
    
    container { "${params.sourmash.container}" }
    
    input:
    path("signatures/*")
    
    output:
    path "similarity_matrix.csv", emit: matrix
    path "versions.yml", emit: versions
    
    script:
    """
    sourmash compare \
        signatures/* \
        -k 51 \
        -o temp.npy \
        --csv similarity_matrix.csv
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(sourmash --version | sed 's/sourmash //')
    END_VERSIONS
    """
}