#!/usr/bin/env python3
"""
Script to update cluster_data_updated.json with optimized coassembly results
from multiple cluster optimization JSON files.

Usage:
    python update_cluster_optimized.py cluster_data.json cluster_144_optimized.json cluster_292_optimized.json
    python update_cluster_optimized.py cluster_data.json -o output.json cluster_*_optimized.json
    python update_cluster_optimized.py cluster_data.json --output new_data.json cluster_*_optimized.json
"""

import json
import sys
import glob
import os
import argparse
from pathlib import Path


def load_json_file(filepath):
    """Load JSON file and return data"""
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            return json.load(f)
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None


def extract_cluster_id_from_filename(filename):
    """Extract cluster ID from filename like 'cluster_144_optimized.json'"""
    basename = os.path.basename(filename)
    if basename.startswith('cluster_') and basename.endswith('_optimized.json'):
        # Extract number between 'cluster_' and '_optimized.json'
        cluster_id = basename[8:-15]  # Remove 'cluster_' and '_optimized.json'
        return cluster_id
    return None


def find_cluster_in_data(cluster_data, cluster_id):
    """Find cluster by ID in the main data structure"""
    for i, cluster in enumerate(cluster_data.get('clusters', [])):
        if cluster.get('cluster_id') == cluster_id:
            return i, cluster
    return None, None


def update_cluster_with_optimized_data(cluster, optimized_data):
    """Update cluster with optimized coassembly results"""
    if 'co_assembly_results' not in cluster:
        cluster['co_assembly_results'] = {}
    
    # Add optimized results
    cluster['co_assembly_results']['optimized'] = {
        "iteration": optimized_data.get("iteration"),
        "test_id": optimized_data.get("test_id"),
        "sag_count": optimized_data.get("sag_count"),
        "sags": optimized_data.get("sags", []),
        "completeness": optimized_data.get("completeness"),
        "contamination": optimized_data.get("contamination"),
        "genome_size": optimized_data.get("genome_size"),
        "contig_n50": optimized_data.get("contig_n50"),
        "contigs_file": optimized_data.get("contigs_file"),
        "output_prefix": optimized_data.get("output_prefix"),
        "result_type": optimized_data.get("result_type"),
        "optimization_summary": optimized_data.get("optimization_summary", {})
    }
    
    return True


def create_backup(filepath):
    """Create a backup of the file"""
    from datetime import datetime
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_path = f"{filepath}.backup_{timestamp}"
    try:
        import shutil
        shutil.copy2(filepath, backup_path)
        print(f"✓ Backup created: {backup_path}")
        return backup_path
    except Exception as e:
        print(f"Warning: Could not create backup: {e}")
        return None


def main():
    parser = argparse.ArgumentParser(
        description='Update cluster data JSON with optimized coassembly results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Output to new file (recommended)
  python update_cluster_optimized.py cluster_data.json -o updated_data.json cluster_*_optimized.json
  
  # Update original file in-place (creates backup)
  python update_cluster_optimized.py cluster_data.json cluster_*_optimized.json
  
  # Specify multiple optimization files
  python update_cluster_optimized.py data.json -o new.json cluster_144_optimized.json cluster_292_optimized.json

Features:
  ✓ Safe output to new file (no modification of original)
  ✓ Automatic backup when updating in-place
  ✓ Cluster ID extraction from filename (cluster_XXX_optimized.json)
  ✓ Wildcard pattern support
  ✓ Progress reporting and error handling
        """
    )
    
    parser.add_argument('main_data_file',
                       help='Path to the main cluster data JSON file')
    parser.add_argument('optimization_files', nargs='+',
                       help='One or more optimization result files (supports wildcards)')
    parser.add_argument('-o', '--output',
                       help='Output file path (if not specified, updates input file in-place with backup)')
    parser.add_argument('--no-backup', action='store_true',
                       help='Skip backup creation when updating in-place (not recommended)')
    
    args = parser.parse_args()
    
    # Get main data file
    main_data_file = args.main_data_file
    
    # Determine output file
    if args.output:
        output_file = args.output
        update_inplace = False
        print(f"Mode: Creating new file")
        print(f"Input:  {main_data_file}")
        print(f"Output: {output_file}")
    else:
        output_file = main_data_file
        update_inplace = True
        print(f"Mode: Updating in-place")
        print(f"File: {main_data_file}")
    
    # Expand glob patterns for optimization files
    input_files = []
    for pattern in args.optimization_files:
        if '*' in pattern or '?' in pattern:
            expanded = glob.glob(pattern)
            if expanded:
                input_files.extend(expanded)
            else:
                print(f"Warning: No files match pattern '{pattern}'")
        else:
            input_files.append(pattern)
    
    if not input_files:
        print("Error: No optimization files found!")
        sys.exit(1)
    
    # Remove duplicates and sort
    input_files = sorted(list(set(input_files)))
    print(f"\nFound {len(input_files)} optimization files to process:")
    for f in input_files:
        print(f"  - {f}")
    print()
    
    if not os.path.exists(main_data_file):
        print(f"Error: Main data file '{main_data_file}' not found!")
        sys.exit(1)
    
    # Create backup if updating in-place
    backup_path = None
    if update_inplace and not args.no_backup:
        backup_path = create_backup(main_data_file)
    
    # Load main data
    print(f"Loading main data from {main_data_file}...")
    main_data = load_json_file(main_data_file)
    if main_data is None:
        sys.exit(1)
    
    # Process each optimized file
    updated_clusters = []
    
    for input_file in input_files:
        print(f"\nProcessing {input_file}...")
        
        # Extract cluster ID from filename
        cluster_id = extract_cluster_id_from_filename(input_file)
        if cluster_id is None:
            print(f"  Warning: Could not extract cluster ID from filename {input_file}")
            continue
        
        print(f"  Extracted cluster ID: {cluster_id}")
        
        # Load optimized data
        optimized_data = load_json_file(input_file)
        if optimized_data is None:
            continue
        
        # Find corresponding cluster in main data
        cluster_index, cluster = find_cluster_in_data(main_data, cluster_id)
        if cluster is None:
            print(f"  Warning: Cluster {cluster_id} not found in main data")
            continue
        
        print(f"  Found cluster {cluster_id} at index {cluster_index}")
        
        # Update cluster with optimized data
        if update_cluster_with_optimized_data(cluster, optimized_data):
            updated_clusters.append(cluster_id)
            print(f"  ✓ Successfully updated cluster {cluster_id}")
            print(f"    - Completeness: {optimized_data.get('completeness', 'N/A')}%")
            print(f"    - Contamination: {optimized_data.get('contamination', 'N/A')}%")
            print(f"    - SAG count: {optimized_data.get('sag_count', 'N/A')}")
        else:
            print(f"  ✗ Failed to update cluster {cluster_id}")
    
    # Save updated data
    if updated_clusters:
        print(f"\nSaving updated data to {output_file}...")
        try:
            # Validate JSON before saving
            json_str = json.dumps(main_data, indent=2, ensure_ascii=False)
            
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(json_str)
            
            print(f"✓ Successfully updated {len(updated_clusters)} clusters: {', '.join(updated_clusters)}")
            print(f"✓ Data saved to: {output_file}")
            
            # Show summary
            print(f"\nSummary:")
            for cluster_id in updated_clusters:
                cluster_index, cluster = find_cluster_in_data(main_data, cluster_id)
                if cluster and 'co_assembly_results' in cluster and 'optimized' in cluster['co_assembly_results']:
                    opt_data = cluster['co_assembly_results']['optimized']
                    print(f"  Cluster {cluster_id}:")
                    print(f"    - Completeness: {opt_data.get('completeness', 'N/A')}%")
                    print(f"    - Contamination: {opt_data.get('contamination', 'N/A')}%")
                    print(f"    - SAGs used: {opt_data.get('sag_count', 'N/A')}")
                    print(f"    - Test ID: {opt_data.get('test_id', 'N/A')}")
            
        except Exception as e:
            print(f"✗ Error saving file: {e}")
            if backup_path:
                print(f"You can restore from backup: {backup_path}")
            sys.exit(1)
    else:
        print("\nNo clusters were updated.")
    
    print(f"\nDone! {len(updated_clusters)} clusters updated successfully.")
    if not update_inplace:
        print(f"Original file '{main_data_file}' remains unchanged.")
        print(f"Updated data saved to '{output_file}'.")


if __name__ == "__main__":
    main()