#!/usr/bin/env python3
"""
Script to update cluster_data.json with GTDB-Tk classification results.

This script reads GTDB-Tk bac120.summary.tsv file and updates the corresponding
clusters in the main data file with taxonomic classification information.

Usage:
    python update_gtdbtk_classification.py <main_data_file> <gtdbtk_summary_file>
    python update_gtdbtk_classification.py <main_data_file> <gtdbtk_summary_file> -o <output_file>
    
Example:
    python update_gtdbtk_classification.py cluster_data.json gtdbtk.bac120.summary.tsv -o updated_data.json
    python update_gtdbtk_classification.py cluster_data.json gtdbtk.bac120.summary.tsv
"""

import json
import sys
import os
import csv
import re
import argparse
from datetime import datetime
from pathlib import Path


def load_json_file(filepath):
    """Load JSON file and return data"""
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            return json.load(f)
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None


def load_gtdbtk_summary(filepath):
    """Load GTDB-Tk summary TSV file and return classification data"""
    classifications = {}
    
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            
            for row in reader:
                user_genome = row.get('user_genome', '')
                classification = row.get('classification', '')
                
                # Extract cluster info from user_genome name
                # Expected format: cluster_1040_CoSAG or similar
                cluster_match = re.search(r'cluster_(\d+)', user_genome)
                if cluster_match:
                    cluster_id = cluster_match.group(1)
                    
                    # Check if this is an optimized result (has iteration info or best_contigs)
                    is_optimized = ('iter_' in user_genome or 
                                  'iteration' in user_genome.lower() or 
                                  'best_contigs' in user_genome)
                    
                    classifications[cluster_id] = {
                        'user_genome': user_genome,
                        'classification': classification,
                        'is_optimized': is_optimized,
                        'full_row': dict(row)  # Store complete row data
                    }
                    
                    print(f"Found classification for cluster {cluster_id}: {user_genome}")
                    if is_optimized:
                        print(f"  -> Detected as optimized result")
                else:
                    print(f"Warning: Could not extract cluster ID from '{user_genome}'")
                    
    except Exception as e:
        print(f"Error loading GTDB-Tk summary file {filepath}: {e}")
        return None
    
    return classifications


def parse_gtdbtk_classification(classification_string):
    """Parse GTDB-Tk classification string into structured taxonomy"""
    taxonomy = {
        'domain': '',
        'phylum': '',
        'class': '',
        'order': '',
        'family': '',
        'genus': '',
        'species': ''
    }
    
    if not classification_string:
        return taxonomy
    
    # Split by semicolon and parse each level
    levels = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    parts = classification_string.split(';')
    
    for i, part in enumerate(parts):
        if i < len(levels):
            # Remove prefix (d__, p__, c__, etc.) and clean up
            clean_part = re.sub(r'^[a-z]__', '', part.strip())
            taxonomy[levels[i]] = clean_part
    
    return taxonomy


def find_cluster_in_data(cluster_data, cluster_id):
    """Find cluster by ID in the main data structure"""
    for i, cluster in enumerate(cluster_data.get('clusters', [])):
        if cluster.get('cluster_id') == cluster_id:
            return i, cluster
    return None, None


def update_cluster_with_gtdbtk_data(cluster, gtdbtk_data, is_optimized=False):
    """Update cluster with GTDB-Tk classification data"""
    
    # Parse the classification string
    taxonomy = parse_gtdbtk_classification(gtdbtk_data['classification'])
    
    # Create classification object
    classification_info = {
        'user_genome': gtdbtk_data['user_genome'],
        'classification_string': gtdbtk_data['classification'],
        'taxonomy': taxonomy,
        'gtdbtk_metadata': {
            'fastani_reference': gtdbtk_data['full_row'].get('fastani_reference', ''),
            'fastani_reference_radius': gtdbtk_data['full_row'].get('fastani_reference_radius', ''),
            'fastani_taxonomy': gtdbtk_data['full_row'].get('fastani_taxonomy', ''),
            'fastani_ani': gtdbtk_data['full_row'].get('fastani_ani', ''),
            'fastani_af': gtdbtk_data['full_row'].get('fastani_af', ''),
            'closest_placement_reference': gtdbtk_data['full_row'].get('closest_placement_reference', ''),
            'closest_placement_radius': gtdbtk_data['full_row'].get('closest_placement_radius', ''),
            'closest_placement_taxonomy': gtdbtk_data['full_row'].get('closest_placement_taxonomy', ''),
            'closest_placement_ani': gtdbtk_data['full_row'].get('closest_placement_ani', ''),
            'closest_placement_af': gtdbtk_data['full_row'].get('closest_placement_af', ''),
            'pplacer_taxonomy': gtdbtk_data['full_row'].get('pplacer_taxonomy', ''),
            'classification_method': gtdbtk_data['full_row'].get('classification_method', ''),
            'note': gtdbtk_data['full_row'].get('note', ''),
            'other_related_references': gtdbtk_data['full_row'].get('other_related_references', ''),
            'msa_percent': gtdbtk_data['full_row'].get('msa_percent', ''),
            'translation_table': gtdbtk_data['full_row'].get('translation_table', ''),
            'red_value': gtdbtk_data['full_row'].get('red_value', ''),
            'warnings': gtdbtk_data['full_row'].get('warnings', '')
        },
        'updated_at': datetime.now().isoformat()
    }
    
    # Ensure co_assembly_results exists
    if 'co_assembly_results' not in cluster:
        cluster['co_assembly_results'] = {}
    
    # Add classification to appropriate section
    if is_optimized and 'optimized' in cluster['co_assembly_results']:
        cluster['co_assembly_results']['optimized']['gtdbtk_classification'] = classification_info
        print(f"  -> Added classification to optimized results")
    else:
        # Add to round_1 or create new section
        if 'round_1' not in cluster['co_assembly_results']:
            cluster['co_assembly_results']['round_1'] = {}
        cluster['co_assembly_results']['round_1']['gtdbtk_classification'] = classification_info
        print(f"  -> Added classification to round_1 results")
    
    return True


def create_backup(filepath):
    """Create a backup of the file"""
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
        description='Update cluster data JSON with GTDB-Tk taxonomic classification results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Output to new file (recommended)
  python update_gtdbtk_classification.py cluster_data.json gtdbtk.bac120.summary.tsv -o updated_data.json
  
  # Update original file in-place (creates backup)
  python update_gtdbtk_classification.py cluster_data.json gtdbtk.bac120.summary.tsv
  
  # Specify custom output location
  python update_gtdbtk_classification.py data.json results/gtdbtk.bac120.summary.tsv --output new_data.json

Features:
  ✓ Safe output to new file (no modification of original)
  ✓ Automatic backup when updating in-place
  ✓ Detects optimized vs regular results
  ✓ Parses full taxonomic hierarchy
  ✓ Preserves all GTDB-Tk metadata
  ✓ Progress reporting and error handling
        """
    )
    
    parser.add_argument('main_data_file',
                       help='Path to the main cluster data JSON file')
    parser.add_argument('gtdbtk_summary_file',
                       help='Path to GTDB-Tk bac120.summary.tsv file')
    parser.add_argument('-o', '--output',
                       help='Output file path (if not specified, updates input file in-place with backup)')
    parser.add_argument('--no-backup', action='store_true',
                       help='Skip backup creation when updating in-place (not recommended)')
    
    args = parser.parse_args()
    
    main_data_file = args.main_data_file
    gtdbtk_summary_file = args.gtdbtk_summary_file
    
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
    
    print(f"GTDB-Tk summary file: {gtdbtk_summary_file}")
    print()
    
    # Check if files exist
    if not os.path.exists(main_data_file):
        print(f"Error: Main data file '{main_data_file}' not found!")
        sys.exit(1)
    
    if not os.path.exists(gtdbtk_summary_file):
        print(f"Error: GTDB-Tk summary file '{gtdbtk_summary_file}' not found!")
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
    
    # Load GTDB-Tk classifications
    print(f"Loading GTDB-Tk classifications from {gtdbtk_summary_file}...")
    classifications = load_gtdbtk_summary(gtdbtk_summary_file)
    if classifications is None:
        sys.exit(1)
    
    print(f"Found {len(classifications)} classifications to process")
    print()
    
    # Process each classification
    updated_clusters = []
    
    for cluster_id, gtdbtk_data in classifications.items():
        print(f"Processing cluster {cluster_id}...")
        
        # Find corresponding cluster in main data
        cluster_index, cluster = find_cluster_in_data(main_data, cluster_id)
        if cluster is None:
            print(f"  Warning: Cluster {cluster_id} not found in main data")
            continue
        
        print(f"  Found cluster {cluster_id} at index {cluster_index}")
        
        # Update cluster with GTDB-Tk data
        if update_cluster_with_gtdbtk_data(cluster, gtdbtk_data, gtdbtk_data['is_optimized']):
            updated_clusters.append(cluster_id)
            print(f"  ✓ Successfully updated cluster {cluster_id}")
            
            # Show taxonomy summary
            taxonomy = parse_gtdbtk_classification(gtdbtk_data['classification'])
            print(f"    - Domain: {taxonomy['domain']}")
            print(f"    - Phylum: {taxonomy['phylum']}")
            print(f"    - Class: {taxonomy['class']}")
            print(f"    - Order: {taxonomy['order']}")
            print(f"    - Family: {taxonomy['family']}")
            print(f"    - Genus: {taxonomy['genus']}")
            print(f"    - Species: {taxonomy['species']}")
        else:
            print(f"  ✗ Failed to update cluster {cluster_id}")
        
        print()
    
    # Save updated data
    if updated_clusters:
        print(f"Saving updated data to {output_file}...")
        try:
            # Validate JSON before saving
            json_str = json.dumps(main_data, indent=2, ensure_ascii=False)
            
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(json_str)
            
            print(f"✓ Successfully updated {len(updated_clusters)} clusters with GTDB-Tk classifications")
            print(f"✓ Data saved to: {output_file}")
            
            # Show summary
            print(f"\nSummary:")
            for cluster_id in updated_clusters:
                cluster_index, cluster = find_cluster_in_data(main_data, cluster_id)
                if cluster and 'co_assembly_results' in cluster:
                    # Check where classification was added
                    has_optimized_class = ('optimized' in cluster['co_assembly_results'] and 
                                         'gtdbtk_classification' in cluster['co_assembly_results']['optimized'])
                    has_round1_class = ('round_1' in cluster['co_assembly_results'] and 
                                      'gtdbtk_classification' in cluster['co_assembly_results']['round_1'])
                    
                    location = "optimized" if has_optimized_class else "round_1" if has_round1_class else "unknown"
                    
                    print(f"  Cluster {cluster_id}: classification added to {location} results")
            
        except Exception as e:
            print(f"✗ Error saving file: {e}")
            if backup_path:
                print(f"You can restore from backup: {backup_path}")
            sys.exit(1)
    else:
        print("No clusters were updated.")
    
    print(f"\nDone! {len(updated_clusters)} clusters updated with GTDB-Tk classifications.")
    if not update_inplace:
        print(f"Original file '{main_data_file}' remains unchanged.")
        print(f"Updated data saved to '{output_file}'.")


if __name__ == "__main__":
    main()