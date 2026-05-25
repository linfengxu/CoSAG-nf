process SOURMASH_MATRIX_CONVERT {
    tag "sourmash_matrix_convert"
    label 'process_medium'

    container 'quay.io/xulf2022/python3.8_bio:v1'

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
