#!/usr/bin/env python3
"""
Create HTML report with embedded data to avoid CORS issues.
This script embeds the JSON data directly into the HTML file.
"""

import json
import os
import sys
import argparse
from datetime import datetime


def load_cluster_data(json_file):
    """Load cluster data from JSON file"""
    try:
        with open(json_file, 'r', encoding='utf-8') as f:
            return json.load(f)
    except Exception as e:
        print(f"Error loading {json_file}: {e}")
        return None


def calculate_summary_stats(data):
    """Calculate summary statistics from cluster data"""
    if not data or 'clusters' not in data:
        return None
    
    clusters = data['clusters']
    
    # Basic counts
    total_sags = sum(cluster.get('cluster_size', 0) for cluster in clusters)
    total_clusters = len(clusters)
    
    # Quality assessment for co-assembly results
    high_quality_cosags = 0
    medium_quality_cosags = 0
    low_quality_cosags = 0
    
    optimized_clusters = 0
    final_quality_genomes = 0
    
    for cluster in clusters:
        # Check for optimization
        co_results = cluster.get('co_assembly_results', {})
        if 'optimized' in co_results:
            optimized_clusters += 1
        
        # Check final genome quality using iteration data if available
        best_result = co_results.get('optimized') or co_results.get('round_1')
        if best_result:
            # Use iteration data if available, otherwise use direct values
            if 'iteration' in best_result:
                comp = best_result.get('completeness', 0)
                cont = best_result.get('contamination', 100)
            else:
                comp = best_result.get('completeness', 0)
                cont = best_result.get('contamination', 100)
            
            # Count quality CoSAGs with separate high and medium categories
            if comp > 90 and cont < 5:
                high_quality_cosags += 1
                final_quality_genomes += 1
            elif comp > 50 and cont < 10:
                medium_quality_cosags += 1
                final_quality_genomes += 1
            else:
                low_quality_cosags += 1
    
    # Cluster size statistics
    cluster_sizes = [cluster.get('cluster_size', 0) for cluster in clusters]
    avg_cluster_size = sum(cluster_sizes) / len(cluster_sizes) if cluster_sizes else 0
    largest_cluster = max(cluster_sizes) if cluster_sizes else 0
    
    return {
        'total_sags': total_sags,
        'high_quality_cosags': high_quality_cosags,
        'medium_quality_cosags': medium_quality_cosags,
        'high_medium_quality_cosags': high_quality_cosags + medium_quality_cosags,
        'low_quality_cosags': low_quality_cosags,
        'total_clusters': total_clusters,
        'optimized_clusters': optimized_clusters,
        'final_quality_genomes': final_quality_genomes,
        'avg_cluster_size': round(avg_cluster_size, 1),
        'largest_cluster': largest_cluster
    }


def generate_sag_table_html(data):
    """Generate SAG table HTML with embedded data"""
    if not data or 'clusters' not in data:
        return "<tr><td colspan='6'>No data available</td></tr>"
    
    rows = []
    for cluster in data['clusters']:
        cluster_id = cluster.get('cluster_id', 'Unknown')
        for member in cluster.get('members', []):
            sag_id = member.get('sag_id', 'Unknown')
            completeness = member.get('completeness', 0)
            contamination = member.get('contamination', 0)
            
            # Determine quality
            if completeness > 90 and contamination < 5:
                quality = 'high'
                quality_text = 'High'
            elif completeness > 50 and contamination < 10:
                quality = 'medium'
                quality_text = 'Medium'
            else:
                quality = 'low'
                quality_text = 'Low'
            
            # Get taxonomy
            taxonomy = 'Unknown'
            co_results = cluster.get('co_assembly_results', {})
            best_result = co_results.get('optimized') or co_results.get('round_1')
            if best_result and 'gtdbtk_classification' in best_result:
                tax = best_result['gtdbtk_classification'].get('taxonomy', {})
                genus = tax.get('genus', '')
                species = tax.get('species', '')
                if genus and species:
                    taxonomy = f"{genus} {species}"
                elif genus:
                    taxonomy = genus
                elif tax.get('family'):
                    taxonomy = tax.get('family')
            
            row = f"""
                <tr>
                    <td>{sag_id}</td>
                    <td>{completeness:.1f}</td>
                    <td>{contamination:.1f}</td>
                    <td>{cluster_id}</td>
                    <td><span class="quality-badge quality-{quality}">{quality_text}</span></td>
                    <td style="font-style: italic;">{taxonomy}</td>
                </tr>
            """
            rows.append(row)
    
    return ''.join(rows)


def generate_cluster_table_html(data):
    """Generate cluster table HTML with embedded data"""
    if not data or 'clusters' not in data:
        return "<tr><td colspan='5'>No data available</td></tr>"
    
    rows = []
    for cluster in data['clusters']:
        cluster_id = cluster.get('cluster_id', 'Unknown')
        cluster_size = cluster.get('cluster_size', 0)
        
        # Get member list
        members = cluster.get('members', [])
        member_list = ', '.join([m.get('sag_id', 'Unknown') for m in members])
        if len(member_list) > 50:
            member_list = member_list[:50] + '...'
        
        # Check optimization status
        co_results = cluster.get('co_assembly_results', {})
        has_optimized = 'optimized' in co_results
        status = 'Optimized' if has_optimized else 'Round 1 Only'
        
        # Get best quality - use iteration data if available
        best_result = co_results.get('optimized') or co_results.get('round_1')
        if best_result:
            if 'iteration' in best_result:
                comp = best_result.get('completeness', 0)
                cont = best_result.get('contamination', 0)
            else:
                comp = best_result.get('completeness', 0)
                cont = best_result.get('contamination', 0)
            quality = f"{comp:.1f}% / {cont:.1f}%"
            
            # Only include if meets quality criteria (completeness>50 and contamination<10)
            if comp > 50 and cont < 10:
                row = f"""
                    <tr>
                        <td>{cluster_id}</td>
                        <td>{cluster_size}</td>
                        <td title="{', '.join([m.get('sag_id', 'Unknown') for m in members])}">{member_list}</td>
                        <td>{status}</td>
                        <td>{quality}</td>
                    </tr>
                """
                rows.append(row)
        else:
            quality = 'N/A'
    
    return ''.join(rows)



def create_embedded_html_report(data, title="SAG Analysis Report", output_file="embedded_report.html"):
    """Create HTML report with embedded data"""
    
    stats = calculate_summary_stats(data)
    current_date = datetime.now().strftime('%Y-%m-%d')
    
    # Generate table content
    sag_table_content = generate_sag_table_html(data)
    cluster_table_content = generate_cluster_table_html(data)
    
    # Embed JSON data as JavaScript variable
    json_data = json.dumps(data, indent=2, ensure_ascii=False)
    
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/plotly.js-dist@2.26.0/plotly.min.js"></script>
    <style>
        /* Embedded CSS */
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}

        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, sans-serif;
            line-height: 1.6;
            color: #333;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
        }}

        .container {{
            max-width: 1400px;
            margin: 0 auto;
            background: white;
            min-height: 100vh;
            box-shadow: 0 0 20px rgba(0,0,0,0.1);
        }}

        /* Header */
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 2rem 0;
            text-align: center;
        }}

        .header-content h1 {{
            font-size: 2.5rem;
            font-weight: 300;
            margin-bottom: 0.5rem;
        }}

        .subtitle {{
            font-size: 1.1rem;
            opacity: 0.9;
            margin-bottom: 1rem;
        }}

        .report-info {{
            font-size: 0.9rem;
            opacity: 0.8;
        }}

        /* Navigation */
        .nav-tabs {{
            background: #f8f9fa;
            border-bottom: 1px solid #dee2e6;
            padding: 0 2rem;
            display: flex;
            justify-content: center;
            overflow-x: auto;
        }}

        .tab-button {{
            background: none;
            border: none;
            padding: 1rem 1.5rem;
            cursor: pointer;
            font-size: 0.9rem;
            color: #6c757d;
            border-bottom: 3px solid transparent;
            transition: all 0.3s ease;
            white-space: nowrap;
        }}

        .tab-button:hover {{
            color: #495057;
            background: rgba(0,0,0,0.05);
        }}

        .tab-button.active {{
            color: #667eea;
            border-bottom-color: #667eea;
            background: white;
        }}

        /* Main Content */
        .main-content {{
            padding: 2rem;
        }}

        .tab-content {{
            display: none;
        }}

        .tab-content.active {{
            display: block;
        }}

        /* Sections */
        section {{
            margin-bottom: 3rem;
        }}

        section h2 {{
            font-size: 1.8rem;
            color: #495057;
            margin-bottom: 1.5rem;
            border-bottom: 2px solid #e9ecef;
            padding-bottom: 0.5rem;
        }}

        /* Stats Grid */
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 1.5rem;
            margin-bottom: 2rem;
        }}

        .stat-card {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 2rem;
            border-radius: 8px;
            text-align: center;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }}

        .stat-number {{
            font-size: 2.5rem;
            font-weight: bold;
            margin-bottom: 0.5rem;
        }}

        .stat-label {{
            font-size: 0.9rem;
            opacity: 0.9;
        }}

        /* Tables */
        .table-container {{
            overflow-x: auto;
            border: 1px solid #e9ecef;
            border-radius: 8px;
            margin: 1.5rem 0;
        }}

        .data-table {{
            width: 100%;
            border-collapse: collapse;
            background: white;
        }}

        .data-table th,
        .data-table td {{
            padding: 0.75rem;
            text-align: left;
            border-bottom: 1px solid #e9ecef;
        }}

        .data-table th {{
            background: #f8f9fa;
            font-weight: 600;
            color: #495057;
            position: sticky;
            top: 0;
            z-index: 10;
        }}

        .data-table tbody tr:hover {{
            background: #f8f9fa;
        }}

        /* Quality badges */
        .quality-badge {{
            padding: 0.25rem 0.5rem;
            border-radius: 12px;
            font-size: 0.8rem;
            font-weight: 500;
            text-transform: uppercase;
        }}

        .quality-high {{
            background: #d4edda;
            color: #155724;
        }}

        .quality-medium {{
            background: #fff3cd;
            color: #856404;
        }}

        .quality-low {{
            background: #f8d7da;
            color: #721c24;
        }}

        /* Improvement indicators */
        .improvement-positive {{
            color: #28a745;
            font-weight: bold;
        }}

        .improvement-negative {{
            color: #dc3545;
            font-weight: bold;
        }}

        .improvement-neutral {{
            color: #6c757d;
        }}

        /* Charts */
        .chart-container {{
            background: white;
            border: 1px solid #e9ecef;
            border-radius: 8px;
            padding: 1.5rem;
            margin: 1.5rem 0;
            box-shadow: 0 2px 4px rgba(0,0,0,0.05);
            height: 400px;
        }}

        /* Controls */
        .controls {{
            display: flex;
            gap: 1rem;
            margin-bottom: 1.5rem;
            flex-wrap: wrap;
        }}

        .controls input,
        .controls select {{
            padding: 0.5rem;
            border: 1px solid #ced4da;
            border-radius: 4px;
            font-size: 0.9rem;
        }}

        .controls input {{
            flex: 1;
            min-width: 200px;
        }}

        /* Responsive design */
        @media (max-width: 768px) {{
            .container {{
                margin: 0;
            }}
            
            .main-content {{
                padding: 1rem;
            }}
            
            .header-content h1 {{
                font-size: 2rem;
            }}
            
            .stats-grid {{
                grid-template-columns: 1fr;
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <!-- Header -->
        <header class="header">
            <div class="header-content">
                <h1>Single-cell Amplified Genome Analysis Report</h1>
                <p class="subtitle">Single-cell Amplified Genome Analysis & Quality Assessment</p>
                <div class="report-info">
                    <span id="report-title">{title}</span> | 
                    <span id="report-date">{current_date}</span> | 
                    <span id="report-version">v1.0.0</span>
                </div>
            </div>
        </header>

        <!-- Navigation -->
        <nav class="nav-tabs">
            <button class="tab-button active" onclick="showTab('overview')">Overview</button>
            <button class="tab-button" onclick="showTab('single-sag')">Single SAG Quality</button>
            <button class="tab-button" onclick="showTab('clustering')">Clustering Analysis</button>
        </nav>

        <!-- Tab Contents -->
        <main class="main-content">
            <!-- Overview Tab -->
            <div id="overview" class="tab-content active">
                <section class="pipeline-overview">
                    <h2>Co-assembly Pipeline Overview</h2>
                    <div class="stats-grid">
                        <div class="stat-card">
                            <div class="stat-number">{stats['total_sags'] if stats else 0}</div>
                            <div class="stat-label">Total SAGs Processed</div>
                        </div>
                        <div class="stat-card">
                            <div class="stat-number">{stats['total_clusters'] if stats else 0}</div>
                            <div class="stat-label">Co-assembly Groups</div>
                        </div>
                        <div class="stat-card">
                            <div class="stat-number">{stats['high_quality_cosags'] if stats else 0}</div>
                            <div class="stat-label">High Quality CoSAGs</div>
                        </div>
                        <div class="stat-card">
                            <div class="stat-number">{stats['medium_quality_cosags'] if stats else 0}</div>
                            <div class="stat-label">Medium Quality CoSAGs</div>
                        </div>
                        <div class="stat-card">
                            <div class="stat-number">{stats['final_quality_genomes'] if stats else 0}</div>
                            <div class="stat-label">Final Quality Genomes</div>
                        </div>
                    </div>
                </section>

                <section class="quality-summary">
                    <h2>Quality Summary</h2>
                    <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 2rem; margin-bottom: 2rem;">
                        <div class="chart-container">
                            <h3 style="text-align: center; margin-bottom: 0.5rem; color: #495057;">Quality Distribution</h3>
                            <p style="text-align: center; margin-bottom: 1rem; color: #6c757d; font-size: 0.9rem; font-style: italic;">
                                Distribution of CoSAGs by quality levels: High (>90% completeness, <5% contamination), 
                                Medium (>50% completeness, <10% contamination), and Low quality
                            </p>
                            <canvas id="quality-overview-chart"></canvas>
                        </div>
                        <div class="chart-container">
                            <h3 style="text-align: center; margin-bottom: 0.5rem; color: #495057;">Phylum Distribution</h3>
                            <p style="text-align: center; margin-bottom: 1rem; color: #6c757d; font-size: 0.9rem; font-style: italic;">
                                Taxonomic diversity at phylum level for high and medium quality CoSAGs, 
                                showing the most abundant bacterial phyla in the dataset
                            </p>
                            <canvas id="phylum-chart"></canvas>
                        </div>
                    </div>
                    
                    <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 2rem; margin-top: 2rem;">
                        <div class="chart-container">
                            <h3 style="text-align: center; margin-bottom: 0.5rem; color: #495057;">Completeness Distribution</h3>
                            <p style="text-align: center; margin-bottom: 1rem; color: #6c757d; font-size: 0.9rem; font-style: italic;">
                                Histogram showing completeness levels of high and medium quality CoSAGs, 
                                indicating genome assembly quality and gene content coverage
                            </p>
                            <canvas id="completeness-chart"></canvas>
                        </div>
                        <div class="chart-container">
                            <h3 style="text-align: center; margin-bottom: 0.5rem; color: #495057;">Contamination Distribution</h3>
                            <p style="text-align: center; margin-bottom: 1rem; color: #6c757d; font-size: 0.9rem; font-style: italic;">
                                Histogram showing contamination levels of high and medium quality CoSAGs, 
                                representing the presence of foreign genetic material
                            </p>
                            <canvas id="contamination-chart"></canvas>
                        </div>
                    </div>
                    
                    <div class="chart-container" style="margin-top: 2rem;">
                        <h3 style="text-align: center; margin-bottom: 0.5rem; color: #495057;">Completeness vs Contamination</h3>
                        <p style="text-align: center; margin-bottom: 1rem; color: #6c757d; font-size: 0.9rem; font-style: italic;">
                            Scatter plot showing the relationship between genome completeness and contamination levels. 
                            Green dots represent high quality CoSAGs, yellow dots represent medium quality CoSAGs
                        </p>
                        <canvas id="scatter-chart"></canvas>
                    </div>
                    
                    <div style="margin-top: 2rem;">
                        <h3 style="margin-bottom: 0.5rem; color: #495057;">High and Medium Quality CoSAGs Details</h3>
                        <p style="margin-bottom: 1rem; color: #6c757d; font-size: 0.9rem; font-style: italic;">
                            Detailed information for all CoSAGs meeting quality criteria (completeness >50%, contamination <10%). 
                            Use search and filters to explore specific clusters or quality levels.
                        </p>
                        <div class="controls">
                            <input type="text" id="cosag-search" placeholder="Search CoSAGs..." onkeyup="filterCoSAGs()">
                            <select id="cosag-sort" onchange="sortCoSAGs()">
                                <option value="cluster_id">Sort by Cluster ID</option>
                                <option value="completeness">Sort by Completeness</option>
                                <option value="contamination">Sort by Contamination</option>
                                <option value="quality">Sort by Quality</option>
                                <option value="size">Sort by Cluster Size</option>
                                <option value="genome_size">Sort by Genome Size</option>
                            </select>
                            <button onclick="downloadCoSAGsCSV()" style="padding: 0.5rem 1rem; background: #667eea; color: white; border: none; border-radius: 4px; cursor: pointer;">Download CSV</button>
                        </div>
                        <div class="table-container">
                            <table id="cosag-table" class="data-table">
                                <thead>
                                    <tr>
                                        <th>Cluster ID</th>
                                        <th>Cluster Size</th>
                                        <th>Completeness (%)</th>
                                        <th>Contamination (%)</th>
                                        <th>Quality</th>
                                        <th>Genome Size (bp)</th>
                                        <th>Taxonomy</th>
                                        <th>Members of SAG ID</th>
                                    </tr>
                                </thead>
                                <tbody id="cosag-table-body">
                                    <!-- Will be populated by JavaScript -->
                                </tbody>
                            </table>
                        </div>
                        <div id="cosag-pagination" style="margin-top: 1rem; text-align: center;">
                            <!-- Pagination controls will be added by JavaScript -->
                        </div>
                    </div>
                </section>
            </div>

            <!-- Single SAG Quality Tab -->
            <div id="single-sag" class="tab-content">
                <section class="sag-quality">
                    <h2>Single SAG Quality Assessment</h2>
                    
                    <!-- Quality Overview -->
                    <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 2rem; margin-bottom: 2rem;">
                        <div class="chart-container">
                            <h3 style="text-align: center; margin-bottom: 0.5rem; color: #495057;">SAG Quality Distribution</h3>
                            <p style="text-align: center; margin-bottom: 1rem; color: #6c757d; font-size: 0.9rem; font-style: italic;">
                                Distribution of individual SAGs by quality levels before co-assembly, 
                                showing the raw quality of single-cell amplified genomes
                            </p>
                            <canvas id="sag-quality-pie-chart"></canvas>
                        </div>
                        <div class="chart-container">
                            <h3 style="text-align: center; margin-bottom: 0.5rem; color: #495057;">SAG Taxonomy Distribution</h3>
                            <p style="text-align: center; margin-bottom: 1rem; color: #6c757d; font-size: 0.9rem; font-style: italic;">
                                Taxonomic classification of individual SAGs based on co-assembly results, 
                                showing the diversity of bacterial taxa in single-cell data
                            </p>
                            <canvas id="sag-taxonomy-chart"></canvas>
                        </div>
                    </div>
                    
                    <!-- Quality Distribution Charts -->
                    <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 2rem; margin-bottom: 2rem;">
                        <div class="chart-container">
                            <h3 style="text-align: center; margin-bottom: 0.5rem; color: #495057;">SAG Completeness Distribution</h3>
                            <p style="text-align: center; margin-bottom: 1rem; color: #6c757d; font-size: 0.9rem; font-style: italic;">
                                Histogram of completeness levels for all individual SAGs, 
                                showing the range of genome assembly quality before co-assembly
                            </p>
                            <canvas id="sag-completeness-chart"></canvas>
                        </div>
                        <div class="chart-container">
                            <h3 style="text-align: center; margin-bottom: 0.5rem; color: #495057;">SAG Contamination Distribution</h3>
                            <p style="text-align: center; margin-bottom: 1rem; color: #6c757d; font-size: 0.9rem; font-style: italic;">
                                Histogram of contamination levels for all individual SAGs, 
                                indicating the presence of foreign genetic material in single-cell data
                            </p>
                            <canvas id="sag-contamination-chart"></canvas>
                        </div>
                    </div>
                    
                    <!-- Scatter Plot -->
                    <div class="chart-container" style="margin-bottom: 2rem;">
                        <h3 style="text-align: center; margin-bottom: 0.5rem; color: #495057;">SAG Quality Scatter Plot</h3>
                        <p style="text-align: center; margin-bottom: 1rem; color: #6c757d; font-size: 0.9rem; font-style: italic;">
                            Scatter plot showing individual SAG quality distribution with color-coded quality levels. 
                            Each point represents one SAG with its completeness and contamination values
                        </p>
                        <canvas id="sag-scatter-chart"></canvas>
                    </div>
                    
                    <!-- Quality vs Cluster Size Analysis -->
                    <div class="chart-container" style="margin-bottom: 2rem;">
                        <h3 style="text-align: center; margin-bottom: 0.5rem; color: #495057;">Average SAG Quality by Cluster Size</h3>
                        <p style="text-align: center; margin-bottom: 1rem; color: #6c757d; font-size: 0.9rem; font-style: italic;">
                            Analysis of how cluster size affects average SAG quality, showing the relationship between 
                            clustering success and individual genome quality metrics
                        </p>
                        <canvas id="quality-vs-cluster-size-chart"></canvas>
                    </div>
                    
                    <!-- SAG Details Table -->
                    <div>
                        <h3 style="margin-bottom: 0.5rem; color: #495057;">Single SAG Details</h3>
                        <p style="margin-bottom: 1rem; color: #6c757d; font-size: 0.9rem; font-style: italic;">
                            Comprehensive table of all individual SAGs with quality metrics and taxonomic information. 
                            Filter by quality level or taxonomy annotation status to explore specific subsets.
                        </p>
                        <div class="controls">
                            <input type="text" id="sag-search" placeholder="Search SAGs..." onkeyup="filterSAGs()">
                            <select id="quality-filter" onchange="filterSAGs()">
                                <option value="all">All SAGs</option>
                                <option value="high">High Quality</option>
                                <option value="medium">Medium Quality</option>
                                <option value="low">Low Quality</option>
                            </select>
                            <select id="taxonomy-filter" onchange="filterSAGs()">
                                <option value="all">All Taxonomy</option>
                                <option value="annotated">With Taxonomy</option>
                                <option value="unknown">Unknown Taxonomy</option>
                            </select>
                            <select id="sag-sort" onchange="sortSAGs()">
                                <option value="sag_id">Sort by SAG ID</option>
                                <option value="completeness">Sort by Completeness</option>
                                <option value="contamination">Sort by Contamination</option>
                                <option value="quality">Sort by Quality</option>
                                <option value="cluster_id">Sort by Cluster ID</option>
                            </select>
                            <button onclick="downloadSAGsCSV()" style="padding: 0.5rem 1rem; background: #667eea; color: white; border: none; border-radius: 4px; cursor: pointer;">Download CSV</button>
                        </div>
                        <div class="table-container">
                            <table id="sag-table" class="data-table">
                                <thead>
                                    <tr>
                                        <th>SAG ID</th>
                                        <th>Completeness (%)</th>
                                        <th>Contamination (%)</th>
                                        <th>Cluster ID</th>
                                        <th>Quality</th>
                                        <th>Taxonomy</th>
                                    </tr>
                                </thead>
                                <tbody id="sag-table-body">
                                    <!-- Will be populated by JavaScript -->
                                </tbody>
                            </table>
                        </div>
                        <div id="sag-pagination" style="margin-top: 1rem; text-align: center;">
                            <!-- Pagination controls will be added by JavaScript -->
                        </div>
                    </div>
                </section>
            </div>

            <!-- Clustering Analysis Tab -->
            <div id="clustering" class="tab-content">
                <section class="clustering-analysis">
                    <h2>Clustering Analysis</h2>
                    
                    <!-- Statistics Overview -->
                    <div class="stats-grid">
                        <div class="stat-card">
                            <div class="stat-number">{stats['total_clusters'] if stats else 0}</div>
                            <div class="stat-label">Total Clusters</div>
                        </div>
                        <div class="stat-card">
                            <div class="stat-number">{stats['avg_cluster_size'] if stats else 0}</div>
                            <div class="stat-label">Average Cluster Size</div>
                        </div>
                        <div class="stat-card">
                            <div class="stat-number">{stats['largest_cluster'] if stats else 0}</div>
                            <div class="stat-label">Largest Cluster</div>
                        </div>
                        <div class="stat-card">
                            <div class="stat-number" id="singleton-clusters">0</div>
                            <div class="stat-label">Singleton Clusters</div>
                        </div>
                    </div>
                    
                    <!-- Clustering Visualization Charts -->
                    <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 2rem; margin: 2rem 0;">
                        <div class="chart-container">
                            <h3 style="text-align: center; margin-bottom: 0.5rem; color: #495057;">Cluster Size Distribution</h3>
                            <p style="text-align: center; margin-bottom: 1rem; color: #6c757d; font-size: 0.9rem; font-style: italic;">
                                Distribution of clusters by size with quality breakdown. Stacked bars show high (green), 
                                medium (yellow), and low (red) quality clusters for each cluster size
                            </p>
                            <canvas id="cluster-size-distribution-chart"></canvas>
                        </div>
                        <div class="chart-container">
                            <h3 style="text-align: center; margin-bottom: 0.5rem; color: #495057;">Clustering Efficiency</h3>
                            <p style="text-align: center; margin-bottom: 1rem; color: #6c757d; font-size: 0.9rem; font-style: italic;">
                                Proportion of SAGs successfully clustered versus remaining as singletons. 
                                Higher clustering efficiency indicates better similarity detection
                            </p>
                            <canvas id="clustering-efficiency-chart"></canvas>
                        </div>
                    </div>
                    
                    <!-- Clustering Impact Analysis -->
                    <div class="chart-container" style="margin: 2rem 0;">
                        <h3 style="text-align: center; margin-bottom: 0.5rem; color: #495057;">Contamination Reduction Through Co-assembly</h3>
                        <p style="text-align: center; margin-bottom: 1rem; color: #6c757d; font-size: 0.9rem; font-style: italic;">
                            Comparison of contamination levels before and after co-assembly optimization, demonstrating the decontamination effect of clustering
                        </p>
                        <canvas id="contamination-comparison-chart"></canvas>
                    </div>
                    
                    <!-- Enhanced Cluster Details Table -->
                    <div style="margin-top: 2rem;">
                        <h3 style="margin-bottom: 0.5rem; color: #495057;">Cluster Details</h3>
                        <p style="margin-bottom: 1rem; color: #6c757d; font-size: 0.9rem; font-style: italic;">
                            Detailed information for each cluster showing co-assembly results and quality improvements. 
                            Quality improvement indicates changes in completeness and contamination after optimization.
                        </p>
                        <div class="controls">
                            <input type="text" id="cluster-search" placeholder="Search clusters..." onkeyup="filterClusters()">
                            <select id="cluster-size-filter" onchange="filterClusters()">
                                <option value="all">All Sizes</option>
                                <option value="1">Singletons (1 SAG)</option>
                                <option value="2">Small (2 SAGs)</option>
                                <option value="3-5">Medium (3-5 SAGs)</option>
                                <option value="6+">Large (6+ SAGs)</option>
                            </select>
                            <select id="cluster-taxonomy-filter" onchange="filterClusters()">
                                <option value="all">All Taxonomy</option>
                                <option value="annotated">With Taxonomy</option>
                                <option value="unknown">Unknown Taxonomy</option>
                            </select>
                            <select id="cluster-sort" onchange="sortClusters()">
                                <option value="cluster_id">Sort by Cluster ID</option>
                                <option value="size">Sort by Size</option>
                                <option value="completeness">Sort by Completeness</option>
                                <option value="contamination">Sort by Contamination</option>
                                <option value="improvement">Sort by Improvement</option>
                            </select>
                            <button onclick="downloadClustersCSV()" style="padding: 0.5rem 1rem; background: #667eea; color: white; border: none; border-radius: 4px; cursor: pointer;">Download CSV</button>
                        </div>
                        <div class="table-container">
                            <table id="cluster-table" class="data-table">
                                <thead>
                                    <tr>
                                        <th>Cluster ID</th>
                                        <th>Size</th>
                                        <th>Members</th>
                                        <th>Completeness (%)</th>
                                        <th>Contamination (%)</th>
                                        <th>Quality Improvement</th>
                                        <th>Taxonomy</th>
                                    </tr>
                                </thead>
                                <tbody id="cluster-table-body">
                                    <!-- Will be populated by JavaScript -->
                                </tbody>
                            </table>
                        </div>
                        <div id="cluster-pagination" style="margin-top: 1rem; text-align: center;">
                            <!-- Pagination controls will be added by JavaScript -->
                        </div>
                    </div>
                </section>
            </div>


        </main>
    </div>

    <script>
        // Embedded data
        const clusterData = {json_data};
        
        // Tab switching functionality
        function showTab(tabName) {{
            // Hide all tab contents
            document.querySelectorAll('.tab-content').forEach(content => {{
                content.classList.remove('active');
            }});

            // Remove active class from all tab buttons
            document.querySelectorAll('.tab-button').forEach(button => {{
                button.classList.remove('active');
            }});

            // Show selected tab content
            const tabContent = document.getElementById(tabName);
            if (tabContent) {{
                tabContent.classList.add('active');
            }}

            // Add active class to selected tab button
            const tabButton = document.querySelector(`[onclick="showTab('${{tabName}}')"]`);
            if (tabButton) {{
                tabButton.classList.add('active');
            }}

            // Load tab-specific content
            if (tabName === 'overview') {{
                createQualityChart();
            }} else if (tabName === 'single-sag') {{
                createSAGCharts();
                populateSAGsTable();
            }} else if (tabName === 'clustering') {{
                createClusteringCharts();
                populateClustersTable();
            }}
        }}

        // Filter SAGs function
        function filterSAGs() {{
            const searchTerm = document.getElementById('sag-search')?.value.toLowerCase() || '';
            const qualityFilter = document.getElementById('quality-filter')?.value || 'all';
            
            const tbody = document.getElementById('sag-table-body');
            if (!tbody) return;

            const rows = tbody.querySelectorAll('tr');
            
            rows.forEach(row => {{
                const sagId = row.cells[0].textContent.toLowerCase();
                const clusterId = row.cells[3].textContent.toLowerCase();
                const taxonomy = row.cells[5].textContent.toLowerCase();
                const qualityBadge = row.querySelector('.quality-badge');
                const quality = qualityBadge ? qualityBadge.textContent.trim().toLowerCase() : '';
                
                const matchesSearch = !searchTerm || 
                    sagId.includes(searchTerm) || 
                    clusterId.includes(searchTerm) || 
                    taxonomy.includes(searchTerm);
                
                const matchesQuality = qualityFilter === 'all' || quality === qualityFilter;
                
                row.style.display = matchesSearch && matchesQuality ? '' : 'none';
            }});
        }}

        // Create quality overview chart
        function createQualityChart() {{
            const ctx = document.getElementById('quality-overview-chart');
            if (!ctx) return;

            // Calculate quality distribution for CoSAGs
            let highCount = 0, mediumCount = 0, lowCount = 0;
            
            if (clusterData && clusterData.clusters) {{
                clusterData.clusters.forEach(cluster => {{
                    const coResults = cluster.co_assembly_results || {{}};
                    const bestResult = coResults.optimized || coResults.round_1;
                    
                    if (bestResult) {{
                        // Use iteration data if available
                        const comp = bestResult.completeness || 0;
                        const cont = bestResult.contamination || 0;
                        
                        if (comp > 90 && cont < 5) {{
                            highCount++;
                        }} else if (comp > 50 && cont < 10) {{
                            mediumCount++;
                        }} else {{
                            lowCount++;
                        }}
                    }}
                }});
            }}

            new Chart(ctx, {{
                type: 'doughnut',
                data: {{
                    labels: ['High Quality CoSAGs', 'Medium Quality CoSAGs', 'Low Quality CoSAGs'],
                    datasets: [{{
                        data: [highCount, mediumCount, lowCount],
                        backgroundColor: ['#28a745', '#ffc107', '#dc3545'],
                        borderWidth: 2,
                        borderColor: '#fff'
                    }}]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {{
                        legend: {{
                            position: 'bottom',
                            labels: {{
                                padding: 20,
                                usePointStyle: true
                            }}
                        }},
                        tooltip: {{
                            callbacks: {{
                                label: function(context) {{
                                    const total = context.dataset.data.reduce((a, b) => a + b, 0);
                                    const percentage = Math.round((context.parsed / total) * 100);
                                    return `${{context.label}}: ${{context.parsed}} (${{percentage}}%)`;
                                }}
                            }}
                        }}
                    }}
                }}
            }});
            
            // Create completeness and contamination charts
            createPhylumChart();
            createCompletenessChart();
            createContaminationChart();
            createScatterChart();
            populateCoSAGsTable();
        }}

        // Create phylum distribution chart
        function createPhylumChart() {{
            const ctx = document.getElementById('phylum-chart');
            if (!ctx) return;

            // Collect phylum data for high and medium quality CoSAGs
            const phylumCounts = {{}};
            
            if (clusterData && clusterData.clusters) {{
                clusterData.clusters.forEach(cluster => {{
                    const coResults = cluster.co_assembly_results || {{}};
                    const bestResult = coResults.optimized || coResults.round_1;
                    
                    if (bestResult) {{
                        const comp = bestResult.completeness || 0;
                        const cont = bestResult.contamination || 0;
                        
                        // Only include high and medium quality CoSAGs
                        if (comp > 50 && cont < 10) {{
                            let phylum = 'Unknown';
                            
                            // Extract phylum from taxonomy
                            if (bestResult.gtdbtk_classification && bestResult.gtdbtk_classification.taxonomy) {{
                                const tax = bestResult.gtdbtk_classification.taxonomy;
                                if (tax.phylum) {{
                                    phylum = tax.phylum;
                                    // Clean up phylum name (remove prefixes like 'p__')
                                    phylum = phylum.replace(/^p__/, '');
                                }}
                            }}
                            
                            phylumCounts[phylum] = (phylumCounts[phylum] || 0) + 1;
                        }}
                    }}
                }});
            }}

            // Convert to arrays for Chart.js
            const phylumEntries = Object.entries(phylumCounts);
            phylumEntries.sort((a, b) => b[1] - a[1]); // Sort by count descending
            
            // Limit to top 8 phyla, group others as "Others"
            const maxPhyla = 8;
            let labels = [];
            let data = [];
            let othersCount = 0;
            
            phylumEntries.forEach((entry, index) => {{
                if (index < maxPhyla) {{
                    labels.push(entry[0]);
                    data.push(entry[1]);
                }} else {{
                    othersCount += entry[1];
                }}
            }});
            
            if (othersCount > 0) {{
                labels.push('Others');
                data.push(othersCount);
            }}

            // Generate colors for phyla
            const colors = [
                '#FF6384', '#36A2EB', '#FFCE56', '#4BC0C0',
                '#9966FF', '#FF9F40', '#E7E9ED', '#71B37C',
                '#8E5EA2', '#F7464A', '#46BFBD', '#FDB45C',
                '#949FB1', '#4D5360', '#AC64AD', '#5CB3CC'
            ];

            new Chart(ctx, {{
                type: 'doughnut',
                data: {{
                    labels: labels,
                    datasets: [{{
                        data: data,
                        backgroundColor: colors.slice(0, labels.length),
                        borderWidth: 2,
                        borderColor: '#fff'
                    }}]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {{
                        legend: {{
                            position: 'bottom',
                            labels: {{
                                padding: 15,
                                usePointStyle: true,
                                font: {{
                                    size: 11
                                }},
                                generateLabels: function(chart) {{
                                    const data = chart.data;
                                    if (data.labels.length && data.datasets.length) {{
                                        return data.labels.map((label, i) => {{
                                            const count = data.datasets[0].data[i];
                                            const total = data.datasets[0].data.reduce((a, b) => a + b, 0);
                                            const percentage = Math.round((count / total) * 100);
                                            return {{
                                                text: `${{label}} (${{count}}, ${{percentage}}%)`,
                                                fillStyle: data.datasets[0].backgroundColor[i],
                                                strokeStyle: data.datasets[0].borderColor,
                                                lineWidth: data.datasets[0].borderWidth,
                                                pointStyle: 'circle',
                                                hidden: false,
                                                index: i
                                            }};
                                        }});
                                    }}
                                    return [];
                                }}
                            }}
                        }},
                        tooltip: {{
                            callbacks: {{
                                label: function(context) {{
                                    const total = context.dataset.data.reduce((a, b) => a + b, 0);
                                    const percentage = Math.round((context.parsed / total) * 100);
                                    return `${{context.label}}: ${{context.parsed}} CoSAGs (${{percentage}}%)`;
                                }}
                            }}
                        }}
                    }}
                }}
            }});
        }}

        // Create completeness distribution chart
        function createCompletenessChart() {{
            const ctx = document.getElementById('completeness-chart');
            if (!ctx) return;

            // Collect completeness data for high and medium quality CoSAGs
            const completenessData = [];
            
            if (clusterData && clusterData.clusters) {{
                clusterData.clusters.forEach(cluster => {{
                    const coResults = cluster.co_assembly_results || {{}};
                    const bestResult = coResults.optimized || coResults.round_1;
                    
                    if (bestResult) {{
                        const comp = bestResult.completeness || 0;
                        const cont = bestResult.contamination || 0;
                        
                        // Only include high and medium quality CoSAGs
                        if (comp > 50 && cont < 10) {{
                            completenessData.push(comp);
                        }}
                    }}
                }});
            }}

            // Create bins for completeness distribution
            const bins = [
                {{ label: '50-60%', min: 50, max: 60, count: 0 }},
                {{ label: '60-70%', min: 60, max: 70, count: 0 }},
                {{ label: '70-80%', min: 70, max: 80, count: 0 }},
                {{ label: '80-90%', min: 80, max: 90, count: 0 }},
                {{ label: '90-100%', min: 90, max: 100, count: 0 }}
            ];

            completenessData.forEach(comp => {{
                bins.forEach(bin => {{
                    if (comp >= bin.min && comp < bin.max) {{
                        bin.count++;
                    }}
                }});
            }});

            new Chart(ctx, {{
                type: 'bar',
                data: {{
                    labels: bins.map(bin => bin.label),
                    datasets: [{{
                        label: 'Number of CoSAGs',
                        data: bins.map(bin => bin.count),
                        backgroundColor: '#667eea',
                        borderColor: '#5a6fd8',
                        borderWidth: 1
                    }}]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {{
                        legend: {{
                            display: false
                        }}
                    }},
                    scales: {{
                        y: {{
                            beginAtZero: true,
                            title: {{
                                display: true,
                                text: 'Number of CoSAGs'
                            }}
                        }},
                        x: {{
                            title: {{
                                display: true,
                                text: 'Completeness Range'
                            }}
                        }}
                    }}
                }}
            }});
        }}

        // Create contamination distribution chart
        function createContaminationChart() {{
            const ctx = document.getElementById('contamination-chart');
            if (!ctx) return;

            // Collect contamination data for high and medium quality CoSAGs
            const contaminationData = [];
            
            if (clusterData && clusterData.clusters) {{
                clusterData.clusters.forEach(cluster => {{
                    const coResults = cluster.co_assembly_results || {{}};
                    const bestResult = coResults.optimized || coResults.round_1;
                    
                    if (bestResult) {{
                        const comp = bestResult.completeness || 0;
                        const cont = bestResult.contamination || 0;
                        
                        // Only include high and medium quality CoSAGs
                        if (comp > 50 && cont < 10) {{
                            contaminationData.push(cont);
                        }}
                    }}
                }});
            }}

            // Create bins for contamination distribution
            const bins = [
                {{ label: '0-2%', min: 0, max: 2, count: 0 }},
                {{ label: '2-4%', min: 2, max: 4, count: 0 }},
                {{ label: '4-6%', min: 4, max: 6, count: 0 }},
                {{ label: '6-8%', min: 6, max: 8, count: 0 }},
                {{ label: '8-10%', min: 8, max: 10, count: 0 }}
            ];

            contaminationData.forEach(cont => {{
                bins.forEach(bin => {{
                    if (cont >= bin.min && cont < bin.max) {{
                        bin.count++;
                    }}
                }});
            }});

            new Chart(ctx, {{
                type: 'bar',
                data: {{
                    labels: bins.map(bin => bin.label),
                    datasets: [{{
                        label: 'Number of CoSAGs',
                        data: bins.map(bin => bin.count),
                        backgroundColor: '#764ba2',
                        borderColor: '#6a4190',
                        borderWidth: 1
                    }}]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {{
                        legend: {{
                            display: false
                        }}
                    }},
                    scales: {{
                        y: {{
                            beginAtZero: true,
                            title: {{
                                display: true,
                                text: 'Number of CoSAGs'
                            }}
                        }},
                        x: {{
                            title: {{
                                display: true,
                                text: 'Contamination Range'
                            }}
                        }}
                    }}
                }}
            }});
        }}

        // Create scatter plot for completeness vs contamination
        function createScatterChart() {{
            const ctx = document.getElementById('scatter-chart');
            if (!ctx) return;

            // Collect data for high and medium quality CoSAGs separately
            const highQualityData = [];
            const mediumQualityData = [];
            
            if (clusterData && clusterData.clusters) {{
                clusterData.clusters.forEach(cluster => {{
                    const coResults = cluster.co_assembly_results || {{}};
                    const bestResult = coResults.optimized || coResults.round_1;
                    
                    if (bestResult) {{
                        const comp = bestResult.completeness || 0;
                        const cont = bestResult.contamination || 0;
                        
                        // Categorize by quality
                        if (comp > 90 && cont < 5) {{
                            // High quality
                            highQualityData.push({{
                                x: cont,
                                y: comp,
                                clusterId: cluster.cluster_id || 'Unknown'
                            }});
                        }} else if (comp > 50 && cont < 10) {{
                            // Medium quality
                            mediumQualityData.push({{
                                x: cont,
                                y: comp,
                                clusterId: cluster.cluster_id || 'Unknown'
                            }});
                        }}
                    }}
                }});
            }}

            new Chart(ctx, {{
                type: 'scatter',
                data: {{
                    datasets: [{{
                        label: 'High Quality CoSAGs',
                        data: highQualityData,
                        backgroundColor: '#28a745',
                        borderColor: '#1e7e34',
                        borderWidth: 1,
                        pointRadius: 6,
                        pointHoverRadius: 8
                    }}, {{
                        label: 'Medium Quality CoSAGs',
                        data: mediumQualityData,
                        backgroundColor: '#ffc107',
                        borderColor: '#e0a800',
                        borderWidth: 1,
                        pointRadius: 5,
                        pointHoverRadius: 7
                    }}]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {{
                        legend: {{
                            display: true,
                            position: 'top',
                            labels: {{
                                padding: 20,
                                usePointStyle: true,
                                font: {{
                                    size: 12
                                }}
                            }}
                        }},
                        tooltip: {{
                            callbacks: {{
                                label: function(context) {{
                                    const point = context.raw;
                                    const quality = context.datasetIndex === 0 ? 'High Quality' : 'Medium Quality';
                                    return `${{quality}} - Cluster ${{point.clusterId}}: ${{point.y.toFixed(1)}}% completeness, ${{point.x.toFixed(1)}}% contamination`;
                                }}
                            }}
                        }}
                    }},
                    scales: {{
                        x: {{
                            title: {{
                                display: true,
                                text: 'Contamination (%)'
                            }},
                            min: 0,
                            max: 10
                        }},
                        y: {{
                            title: {{
                                display: true,
                                text: 'Completeness (%)'
                            }},
                            min: 50,
                            max: 100
                        }}
                    }}
                }}
            }});
        }}

        // Global variables for table management
        let allCoSAGsData = [];
        let filteredCoSAGsData = [];
        let currentPage = 1;
        const itemsPerPage = 20;

        // Populate CoSAGs table
        function populateCoSAGsTable() {{
            allCoSAGsData = [];
            
            if (clusterData && clusterData.clusters) {{
                clusterData.clusters.forEach(cluster => {{
                    const coResults = cluster.co_assembly_results || {{}};
                    const bestResult = coResults.optimized || coResults.round_1;
                    
                    if (bestResult) {{
                        const comp = bestResult.completeness || 0;
                        const cont = bestResult.contamination || 0;
                        
                        // Only include high and medium quality CoSAGs
                        if (comp > 50 && cont < 10) {{
                            // Get taxonomy - use full classification string
                            let taxonomy = 'Unknown';
                            if (bestResult.gtdbtk_classification && bestResult.gtdbtk_classification.classification_string) {{
                                taxonomy = bestResult.gtdbtk_classification.classification_string;
                            }}
                            
                            // Get member list
                            const members = cluster.members || [];
                            const memberList = members.map(m => m.sag_id || 'Unknown').join(', ');
                            
                            // Get genome size
                            const genomeSize = bestResult.genome_size || 0;
                            
                            // Determine quality level
                            let quality = 'Medium';
                            if (comp > 90 && cont < 5) {{
                                quality = 'High';
                            }}
                            
                            allCoSAGsData.push({{
                                cluster_id: cluster.cluster_id || 'Unknown',
                                cluster_size: cluster.cluster_size || 0,
                                completeness: comp,
                                contamination: cont,
                                quality: quality,
                                genome_size: genomeSize,
                                taxonomy: taxonomy,
                                members: memberList
                            }});
                        }}
                    }}
                }});
            }}
            
            filteredCoSAGsData = [...allCoSAGsData];
            displayCoSAGsTable();
        }}

        // Display CoSAGs table with pagination
        function displayCoSAGsTable() {{
            const tbody = document.getElementById('cosag-table-body');
            if (!tbody) return;
            
            const startIndex = (currentPage - 1) * itemsPerPage;
            const endIndex = startIndex + itemsPerPage;
            const pageData = filteredCoSAGsData.slice(startIndex, endIndex);
            
            tbody.innerHTML = pageData.map(item => `
                <tr>
                    <td>${{item.cluster_id}}</td>
                    <td>${{item.cluster_size}}</td>
                    <td>${{item.completeness.toFixed(1)}}</td>
                    <td>${{item.contamination.toFixed(1)}}</td>
                    <td><span class="quality-badge quality-${{item.quality.toLowerCase()}}">${{item.quality}}</span></td>
                    <td>${{item.genome_size.toLocaleString()}}</td>
                    <td style="font-style: italic; max-width: 300px; word-wrap: break-word; font-size: 0.85rem;" title="${{item.taxonomy}}">${{item.taxonomy}}</td>
                    <td title="${{item.members}}">${{item.members.length > 50 ? item.members.substring(0, 50) + '...' : item.members}}</td>
                </tr>
            `).join('');
            
            updatePagination();
        }}

        // Update pagination controls
        function updatePagination() {{
            const paginationDiv = document.getElementById('cosag-pagination');
            if (!paginationDiv) return;
            
            const totalPages = Math.ceil(filteredCoSAGsData.length / itemsPerPage);
            const totalItems = filteredCoSAGsData.length;
            const startItem = (currentPage - 1) * itemsPerPage + 1;
            const endItem = Math.min(currentPage * itemsPerPage, totalItems);
            
            let paginationHTML = `
                <div style="display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 1rem;">
                    <div>Showing ${{startItem}}-${{endItem}} of ${{totalItems}} CoSAGs</div>
                    <div style="display: flex; gap: 0.5rem; align-items: center;">
            `;
            
            // Previous button
            if (currentPage > 1) {{
                paginationHTML += `<button onclick="changePage(${{currentPage - 1}})" style="padding: 0.5rem 1rem; background: #667eea; color: white; border: none; border-radius: 4px; cursor: pointer;">Previous</button>`;
            }}
            
            // Page numbers
            const maxVisiblePages = 5;
            let startPage = Math.max(1, currentPage - Math.floor(maxVisiblePages / 2));
            let endPage = Math.min(totalPages, startPage + maxVisiblePages - 1);
            
            if (endPage - startPage + 1 < maxVisiblePages) {{
                startPage = Math.max(1, endPage - maxVisiblePages + 1);
            }}
            
            for (let i = startPage; i <= endPage; i++) {{
                const isActive = i === currentPage;
                paginationHTML += `<button onclick="changePage(${{i}})" style="padding: 0.5rem 1rem; background: ${{isActive ? '#495057' : '#f8f9fa'}}; color: ${{isActive ? 'white' : '#495057'}}; border: 1px solid #dee2e6; border-radius: 4px; cursor: pointer;">${{i}}</button>`;
            }}
            
            // Next button
            if (currentPage < totalPages) {{
                paginationHTML += `<button onclick="changePage(${{currentPage + 1}})" style="padding: 0.5rem 1rem; background: #667eea; color: white; border: none; border-radius: 4px; cursor: pointer;">Next</button>`;
            }}
            
            paginationHTML += `
                    </div>
                </div>
            `;
            
            paginationDiv.innerHTML = paginationHTML;
        }}

        // Change page
        function changePage(page) {{
            currentPage = page;
            displayCoSAGsTable();
        }}

        // Filter CoSAGs
        function filterCoSAGs() {{
            const searchTerm = document.getElementById('cosag-search')?.value.toLowerCase() || '';
            
            filteredCoSAGsData = allCoSAGsData.filter(item => {{
                return item.cluster_id.toLowerCase().includes(searchTerm) ||
                       item.taxonomy.toLowerCase().includes(searchTerm) ||
                       item.members.toLowerCase().includes(searchTerm) ||
                       item.genome_size.toString().includes(searchTerm);
            }});
            
            currentPage = 1;
            displayCoSAGsTable();
        }}

        // Sort CoSAGs
        function sortCoSAGs() {{
            const sortBy = document.getElementById('cosag-sort')?.value || 'cluster_id';
            
            filteredCoSAGsData.sort((a, b) => {{
                switch (sortBy) {{
                    case 'completeness':
                        return b.completeness - a.completeness;
                    case 'contamination':
                        return a.contamination - b.contamination;
                    case 'quality':
                        // High quality first, then Medium
                        const qualityOrder = {{ 'High': 2, 'Medium': 1 }};
                        return qualityOrder[b.quality] - qualityOrder[a.quality];
                    case 'size':
                        return b.cluster_size - a.cluster_size;
                    case 'genome_size':
                        return b.genome_size - a.genome_size;
                    case 'cluster_id':
                    default:
                        return a.cluster_id.localeCompare(b.cluster_id);
                }}
            }});
            
            currentPage = 1;
            displayCoSAGsTable();
        }}

        // Download CSV
        function downloadCoSAGsCSV() {{
            const headers = ['Cluster ID', 'Cluster Size', 'Completeness (%)', 'Contamination (%)', 'Quality', 'Genome Size (bp)', 'Taxonomy', 'Members'];
            const csvContent = [
                headers.join(','),
                ...filteredCoSAGsData.map(item => [
                    item.cluster_id,
                    item.cluster_size,
                    item.completeness.toFixed(1),
                    item.contamination.toFixed(1),
                    item.quality,
                    item.genome_size,
                    `"${{item.taxonomy}}"`,
                    `"${{item.members}}"`
                ].join(','))
            ].join('\\n');
            
            const blob = new Blob([csvContent], {{ type: 'text/csv;charset=utf-8;' }});
            const link = document.createElement('a');
            const url = URL.createObjectURL(blob);
            link.setAttribute('href', url);
            link.setAttribute('download', 'high_medium_quality_cosags.csv');
            link.style.visibility = 'hidden';
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
        }}

        // SAG-related functions
        let allSAGsData = [];
        let filteredSAGsData = [];
        let currentSAGPage = 1;
        const sagItemsPerPage = 25;

        // Create SAG charts
        function createSAGCharts() {{
            createSAGQualityPieChart();
            createSAGTaxonomyChart();
            createSAGCompletenessChart();
            createSAGContaminationChart();
            createSAGScatterChart();
            createQualityVsClusterSizeChart();
        }}

        // Create SAG quality pie chart
        function createSAGQualityPieChart() {{
            const ctx = document.getElementById('sag-quality-pie-chart');
            if (!ctx) return;

            // Calculate quality distribution for all SAGs
            let highCount = 0, mediumCount = 0, lowCount = 0;
            
            if (clusterData && clusterData.clusters) {{
                clusterData.clusters.forEach(cluster => {{
                    cluster.members.forEach(member => {{
                        const comp = member.completeness || 0;
                        const cont = member.contamination || 0;
                        
                        if (comp > 90 && cont < 5) {{
                            highCount++;
                        }} else if (comp > 50 && cont < 10) {{
                            mediumCount++;
                        }} else {{
                            lowCount++;
                        }}
                    }});
                }});
            }}

            new Chart(ctx, {{
                type: 'doughnut',
                data: {{
                    labels: ['High Quality SAGs', 'Medium Quality SAGs', 'Low Quality SAGs'],
                    datasets: [{{
                        data: [highCount, mediumCount, lowCount],
                        backgroundColor: ['#28a745', '#ffc107', '#dc3545'],
                        borderWidth: 2,
                        borderColor: '#fff'
                    }}]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {{
                        legend: {{
                            position: 'bottom',
                            labels: {{
                                padding: 20,
                                usePointStyle: true
                            }}
                        }},
                        tooltip: {{
                            callbacks: {{
                                label: function(context) {{
                                    const total = context.dataset.data.reduce((a, b) => a + b, 0);
                                    const percentage = Math.round((context.parsed / total) * 100);
                                    return `${{context.label}}: ${{context.parsed}} SAGs (${{percentage}}%)`;
                                }}
                            }}
                        }}
                    }}
                }}
            }});
        }}

        // Create SAG taxonomy distribution chart
        function createSAGTaxonomyChart() {{
            const ctx = document.getElementById('sag-taxonomy-chart');
            if (!ctx) return;

            // Collect taxonomy data for all SAGs
            const phylumCounts = {{}};
            
            if (clusterData && clusterData.clusters) {{
                clusterData.clusters.forEach(cluster => {{
                    // Get taxonomy from co-assembly results for this cluster
                    const coResults = cluster.co_assembly_results || {{}};
                    const bestResult = coResults.optimized || coResults.round_1;
                    
                    let phylum = 'Unknown';
                    if (bestResult && bestResult.gtdbtk_classification && bestResult.gtdbtk_classification.taxonomy) {{
                        const tax = bestResult.gtdbtk_classification.taxonomy;
                        if (tax.phylum) {{
                            phylum = tax.phylum;
                            // Clean up phylum name (remove prefixes like 'p__')
                            phylum = phylum.replace(/^p__/, '');
                        }}
                    }}
                    
                    // Count each SAG in this cluster with the cluster's taxonomy
                    cluster.members.forEach(member => {{
                        phylumCounts[phylum] = (phylumCounts[phylum] || 0) + 1;
                    }});
                }});
            }}

            // Convert to arrays for Chart.js
            const phylumEntries = Object.entries(phylumCounts);
            phylumEntries.sort((a, b) => b[1] - a[1]); // Sort by count descending
            
            // Limit to top 8 phyla, group others as "Others"
            const maxPhyla = 8;
            let labels = [];
            let data = [];
            let othersCount = 0;
            
            phylumEntries.forEach((entry, index) => {{
                if (index < maxPhyla) {{
                    labels.push(entry[0]);
                    data.push(entry[1]);
                }} else {{
                    othersCount += entry[1];
                }}
            }});
            
            if (othersCount > 0) {{
                labels.push('Others');
                data.push(othersCount);
            }}

            // Generate colors for phyla
            const colors = [
                '#FF6384', '#36A2EB', '#FFCE56', '#4BC0C0',
                '#9966FF', '#FF9F40', '#E7E9ED', '#71B37C',
                '#8E5EA2', '#F7464A', '#46BFBD', '#FDB45C'
            ];

            new Chart(ctx, {{
                type: 'doughnut',
                data: {{
                    labels: labels,
                    datasets: [{{
                        data: data,
                        backgroundColor: colors.slice(0, labels.length),
                        borderWidth: 2,
                        borderColor: '#fff'
                    }}]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {{
                        legend: {{
                            position: 'bottom',
                            labels: {{
                                padding: 15,
                                usePointStyle: true,
                                font: {{
                                    size: 11
                                }},
                                generateLabels: function(chart) {{
                                    const data = chart.data;
                                    if (data.labels.length && data.datasets.length) {{
                                        return data.labels.map((label, i) => {{
                                            const count = data.datasets[0].data[i];
                                            const total = data.datasets[0].data.reduce((a, b) => a + b, 0);
                                            const percentage = Math.round((count / total) * 100);
                                            return {{
                                                text: `${{label}} (${{count}}, ${{percentage}}%)`,
                                                fillStyle: data.datasets[0].backgroundColor[i],
                                                strokeStyle: data.datasets[0].borderColor,
                                                lineWidth: data.datasets[0].borderWidth,
                                                pointStyle: 'circle',
                                                hidden: false,
                                                index: i
                                            }};
                                        }});
                                    }}
                                    return [];
                                }}
                            }}
                        }},
                        tooltip: {{
                            callbacks: {{
                                label: function(context) {{
                                    const total = context.dataset.data.reduce((a, b) => a + b, 0);
                                    const percentage = Math.round((context.parsed / total) * 100);
                                    return `${{context.label}}: ${{context.parsed}} SAGs (${{percentage}}%)`;
                                }}
                            }}
                        }}
                    }}
                }}
            }});
        }}

        // Create SAG completeness distribution chart
        function createSAGCompletenessChart() {{
            const ctx = document.getElementById('sag-completeness-chart');
            if (!ctx) return;

            // Collect completeness data for all SAGs
            const completenessData = [];
            
            if (clusterData && clusterData.clusters) {{
                clusterData.clusters.forEach(cluster => {{
                    cluster.members.forEach(member => {{
                        const comp = member.completeness || 0;
                        completenessData.push(comp);
                    }});
                }});
            }}

            // Create bins for completeness distribution
            const bins = [
                {{ label: '0-10%', min: 0, max: 10, count: 0 }},
                {{ label: '10-20%', min: 10, max: 20, count: 0 }},
                {{ label: '20-30%', min: 20, max: 30, count: 0 }},
                {{ label: '30-40%', min: 30, max: 40, count: 0 }},
                {{ label: '40-50%', min: 40, max: 50, count: 0 }},
                {{ label: '50-60%', min: 50, max: 60, count: 0 }},
                {{ label: '60-70%', min: 60, max: 70, count: 0 }},
                {{ label: '70-80%', min: 70, max: 80, count: 0 }},
                {{ label: '80-90%', min: 80, max: 90, count: 0 }},
                {{ label: '90-100%', min: 90, max: 100, count: 0 }}
            ];

            completenessData.forEach(comp => {{
                bins.forEach(bin => {{
                    if (comp >= bin.min && comp < bin.max) {{
                        bin.count++;
                    }}
                }});
            }});

            new Chart(ctx, {{
                type: 'bar',
                data: {{
                    labels: bins.map(bin => bin.label),
                    datasets: [{{
                        label: 'Number of SAGs',
                        data: bins.map(bin => bin.count),
                        backgroundColor: '#667eea',
                        borderColor: '#5a6fd8',
                        borderWidth: 1
                    }}]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {{
                        legend: {{
                            display: false
                        }}
                    }},
                    scales: {{
                        y: {{
                            beginAtZero: true,
                            title: {{
                                display: true,
                                text: 'Number of SAGs'
                            }}
                        }},
                        x: {{
                            title: {{
                                display: true,
                                text: 'Completeness Range'
                            }}
                        }}
                    }}
                }}
            }});
        }}

        // Create SAG contamination distribution chart
        function createSAGContaminationChart() {{
            const ctx = document.getElementById('sag-contamination-chart');
            if (!ctx) return;

            // Collect contamination data for all SAGs
            const contaminationData = [];
            
            if (clusterData && clusterData.clusters) {{
                clusterData.clusters.forEach(cluster => {{
                    cluster.members.forEach(member => {{
                        const cont = member.contamination || 0;
                        contaminationData.push(cont);
                    }});
                }});
            }}

            // Create bins for contamination distribution
            const bins = [
                {{ label: '0-2%', min: 0, max: 2, count: 0 }},
                {{ label: '2-5%', min: 2, max: 5, count: 0 }},
                {{ label: '5-10%', min: 5, max: 10, count: 0 }},
                {{ label: '10-15%', min: 10, max: 15, count: 0 }},
                {{ label: '15-20%', min: 15, max: 20, count: 0 }},
                {{ label: '20%+', min: 20, max: 100, count: 0 }}
            ];

            contaminationData.forEach(cont => {{
                bins.forEach(bin => {{
                    if (cont >= bin.min && cont < bin.max) {{
                        bin.count++;
                    }}
                }});
            }});

            new Chart(ctx, {{
                type: 'bar',
                data: {{
                    labels: bins.map(bin => bin.label),
                    datasets: [{{
                        label: 'Number of SAGs',
                        data: bins.map(bin => bin.count),
                        backgroundColor: '#764ba2',
                        borderColor: '#6a4190',
                        borderWidth: 1
                    }}]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {{
                        legend: {{
                            display: false
                        }}
                    }},
                    scales: {{
                        y: {{
                            beginAtZero: true,
                            title: {{
                                display: true,
                                text: 'Number of SAGs'
                            }}
                        }},
                        x: {{
                            title: {{
                                display: true,
                                text: 'Contamination Range'
                            }}
                        }}
                    }}
                }}
            }});
        }}

        // Create SAG scatter plot
        function createSAGScatterChart() {{
            const ctx = document.getElementById('sag-scatter-chart');
            if (!ctx) return;

            // Collect data for all SAGs with quality color coding
            const scatterData = [];
            
            if (clusterData && clusterData.clusters) {{
                clusterData.clusters.forEach(cluster => {{
                    cluster.members.forEach(member => {{
                        const comp = member.completeness || 0;
                        const cont = member.contamination || 0;
                        
                        // Determine quality and color
                        let quality, color;
                        if (comp > 90 && cont < 5) {{
                            quality = 'High';
                            color = '#28a745';
                        }} else if (comp > 50 && cont < 10) {{
                            quality = 'Medium';
                            color = '#ffc107';
                        }} else {{
                            quality = 'Low';
                            color = '#dc3545';
                        }}
                        
                        scatterData.push({{
                            x: cont,
                            y: comp,
                            sagId: member.sag_id || 'Unknown',
                            clusterId: cluster.cluster_id || 'Unknown',
                            quality: quality,
                            backgroundColor: color
                        }});
                    }});
                }});
            }}

            new Chart(ctx, {{
                type: 'scatter',
                data: {{
                    datasets: [{{
                        label: 'SAGs',
                        data: scatterData,
                        backgroundColor: function(context) {{
                            return context.raw.backgroundColor;
                        }},
                        borderColor: function(context) {{
                            return context.raw.backgroundColor;
                        }},
                        borderWidth: 1,
                        pointRadius: 4,
                        pointHoverRadius: 6
                    }}]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {{
                        legend: {{
                            display: false
                        }},
                        tooltip: {{
                            callbacks: {{
                                label: function(context) {{
                                    const point = context.raw;
                                    return `${{point.sagId}} (Cluster ${{point.clusterId}}): ${{point.y.toFixed(1)}}% completeness, ${{point.x.toFixed(1)}}% contamination (${{point.quality}} Quality)`;
                                }}
                            }}
                        }}
                    }},
                    scales: {{
                        x: {{
                            title: {{
                                display: true,
                                text: 'Contamination (%)'
                            }},
                            min: 0
                        }},
                        y: {{
                            title: {{
                                display: true,
                                text: 'Completeness (%)'
                            }},
                            min: 0,
                            max: 100
                        }}
                    }}
                }}
            }});
        }}

        // Create cluster size distribution chart
        function createClusterSizeChart() {{
            const ctx = document.getElementById('cluster-size-chart');
            if (!ctx) return;

            // Collect cluster size data
            const clusterSizes = [];
            
            if (clusterData && clusterData.clusters) {{
                clusterData.clusters.forEach(cluster => {{
                    const size = cluster.cluster_size || 0;
                    clusterSizes.push(size);
                }});
            }}

            // Create bins for cluster size distribution
            const bins = [
                {{ label: '1 SAG', min: 1, max: 1, count: 0 }},
                {{ label: '2 SAGs', min: 2, max: 2, count: 0 }},
                {{ label: '3 SAGs', min: 3, max: 3, count: 0 }},
                {{ label: '4 SAGs', min: 4, max: 4, count: 0 }},
                {{ label: '5 SAGs', min: 5, max: 5, count: 0 }},
                {{ label: '6+ SAGs', min: 6, max: 100, count: 0 }}
            ];

            clusterSizes.forEach(size => {{
                bins.forEach(bin => {{
                    if (size >= bin.min && size <= bin.max) {{
                        bin.count++;
                    }}
                }});
            }});

            new Chart(ctx, {{
                type: 'bar',
                data: {{
                    labels: bins.map(bin => bin.label),
                    datasets: [{{
                        label: 'Number of Clusters',
                        data: bins.map(bin => bin.count),
                        backgroundColor: '#36A2EB',
                        borderColor: '#2E8BC0',
                        borderWidth: 1
                    }}]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {{
                        legend: {{
                            display: false
                        }}
                    }},
                    scales: {{
                        y: {{
                            beginAtZero: true,
                            title: {{
                                display: true,
                                text: 'Number of Clusters'
                            }}
                        }},
                        x: {{
                            title: {{
                                display: true,
                                text: 'Cluster Size'
                            }}
                        }}
                    }}
                }}
            }});
        }}

        // Create quality vs cluster size analysis chart
        function createQualityVsClusterSizeChart() {{
            const ctx = document.getElementById('quality-vs-cluster-size-chart');
            if (!ctx) return;

            // Collect data for quality vs cluster size analysis
            const clusterData_analysis = [];
            
            if (clusterData && clusterData.clusters) {{
                clusterData.clusters.forEach(cluster => {{
                    const size = cluster.cluster_size || 0;
                    let totalCompleteness = 0;
                    let totalContamination = 0;
                    let memberCount = 0;
                    
                    cluster.members.forEach(member => {{
                        totalCompleteness += member.completeness || 0;
                        totalContamination += member.contamination || 0;
                        memberCount++;
                    }});
                    
                    if (memberCount > 0) {{
                        clusterData_analysis.push({{
                            size: size,
                            avgCompleteness: totalCompleteness / memberCount,
                            avgContamination: totalContamination / memberCount
                        }});
                    }}
                }});
            }}

            // Group by cluster size and calculate averages
            const sizeGroups = {{}};
            clusterData_analysis.forEach(item => {{
                const sizeKey = item.size > 5 ? '6+' : item.size.toString();
                if (!sizeGroups[sizeKey]) {{
                    sizeGroups[sizeKey] = {{ completeness: [], contamination: [] }};
                }}
                sizeGroups[sizeKey].completeness.push(item.avgCompleteness);
                sizeGroups[sizeKey].contamination.push(item.avgContamination);
            }});

            const labels = Object.keys(sizeGroups).sort((a, b) => {{
                if (a === '6+') return 1;
                if (b === '6+') return -1;
                return parseInt(a) - parseInt(b);
            }});

            const avgCompleteness = labels.map(label => {{
                const values = sizeGroups[label].completeness;
                return values.reduce((a, b) => a + b, 0) / values.length;
            }});

            const avgContamination = labels.map(label => {{
                const values = sizeGroups[label].contamination;
                return values.reduce((a, b) => a + b, 0) / values.length;
            }});

            new Chart(ctx, {{
                type: 'line',
                data: {{
                    labels: labels.map(l => l + ' SAG' + (l === '1' ? '' : 's')),
                    datasets: [{{
                        label: 'Average Completeness (%)',
                        data: avgCompleteness,
                        borderColor: '#28a745',
                        backgroundColor: 'rgba(40, 167, 69, 0.1)',
                        yAxisID: 'y'
                    }}, {{
                        label: 'Average Contamination (%)',
                        data: avgContamination,
                        borderColor: '#dc3545',
                        backgroundColor: 'rgba(220, 53, 69, 0.1)',
                        yAxisID: 'y1'
                    }}]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    interaction: {{
                        mode: 'index',
                        intersect: false,
                    }},
                    plugins: {{
                        legend: {{
                            position: 'top'
                        }}
                    }},
                    scales: {{
                        x: {{
                            title: {{
                                display: true,
                                text: 'Cluster Size'
                            }}
                        }},
                        y: {{
                            type: 'linear',
                            display: true,
                            position: 'left',
                            title: {{
                                display: true,
                                text: 'Average Completeness (%)'
                            }},
                            min: 0
                        }},
                        y1: {{
                            type: 'linear',
                            display: true,
                            position: 'right',
                            title: {{
                                display: true,
                                text: 'Average Contamination (%)'
                            }},
                            min: 0,
                            grid: {{
                                drawOnChartArea: false,
                            }},
                        }}
                    }}
                }}
            }});
        }}

        // Populate SAGs table
        function populateSAGsTable() {{
            allSAGsData = [];
            
            if (clusterData && clusterData.clusters) {{
                clusterData.clusters.forEach(cluster => {{
                    // Get taxonomy from co-assembly results - use full classification string
                    const coResults = cluster.co_assembly_results || {{}};
                    const bestResult = coResults.optimized || coResults.round_1;
                    let taxonomy = 'Unknown';
                    
                    if (bestResult && bestResult.gtdbtk_classification && bestResult.gtdbtk_classification.classification_string) {{
                        taxonomy = bestResult.gtdbtk_classification.classification_string;
                    }}
                    
                    cluster.members.forEach(member => {{
                        const comp = member.completeness || 0;
                        const cont = member.contamination || 0;
                        
                        // Determine quality
                        let quality;
                        if (comp > 90 && cont < 5) {{
                            quality = 'High';
                        }} else if (comp > 50 && cont < 10) {{
                            quality = 'Medium';
                        }} else {{
                            quality = 'Low';
                        }}
                        
                        allSAGsData.push({{
                            sag_id: member.sag_id || 'Unknown',
                            completeness: comp,
                            contamination: cont,
                            cluster_id: cluster.cluster_id || 'Unknown',
                            quality: quality,
                            taxonomy: taxonomy
                        }});
                    }});
                }});
            }}
            
            filteredSAGsData = [...allSAGsData];
            displaySAGsTable();
        }}

        // Display SAGs table with pagination
        function displaySAGsTable() {{
            const tbody = document.getElementById('sag-table-body');
            if (!tbody) return;
            
            const startIndex = (currentSAGPage - 1) * sagItemsPerPage;
            const endIndex = startIndex + sagItemsPerPage;
            const pageData = filteredSAGsData.slice(startIndex, endIndex);
            
            tbody.innerHTML = pageData.map(item => `
                <tr>
                    <td>${{item.sag_id}}</td>
                    <td>${{item.completeness.toFixed(1)}}</td>
                    <td>${{item.contamination.toFixed(1)}}</td>
                    <td>${{item.cluster_id}}</td>
                    <td><span class="quality-badge quality-${{item.quality.toLowerCase()}}">${{item.quality}}</span></td>
                    <td style="font-style: italic; max-width: 300px; word-wrap: break-word; font-size: 0.85rem;" title="${{item.taxonomy}}">${{item.taxonomy}}</td>
                </tr>
            `).join('');
            
            updateSAGPagination();
        }}

        // Update SAG pagination controls
        function updateSAGPagination() {{
            const paginationDiv = document.getElementById('sag-pagination');
            if (!paginationDiv) return;
            
            const totalPages = Math.ceil(filteredSAGsData.length / sagItemsPerPage);
            const totalItems = filteredSAGsData.length;
            const startItem = (currentSAGPage - 1) * sagItemsPerPage + 1;
            const endItem = Math.min(currentSAGPage * sagItemsPerPage, totalItems);
            
            let paginationHTML = `
                <div style="display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 1rem;">
                    <div>Showing ${{startItem}}-${{endItem}} of ${{totalItems}} SAGs</div>
                    <div style="display: flex; gap: 0.5rem; align-items: center;">
            `;
            
            // Previous button
            if (currentSAGPage > 1) {{
                paginationHTML += `<button onclick="changeSAGPage(${{currentSAGPage - 1}})" style="padding: 0.5rem 1rem; background: #667eea; color: white; border: none; border-radius: 4px; cursor: pointer;">Previous</button>`;
            }}
            
            // Page numbers
            const maxVisiblePages = 5;
            let startPage = Math.max(1, currentSAGPage - Math.floor(maxVisiblePages / 2));
            let endPage = Math.min(totalPages, startPage + maxVisiblePages - 1);
            
            if (endPage - startPage + 1 < maxVisiblePages) {{
                startPage = Math.max(1, endPage - maxVisiblePages + 1);
            }}
            
            for (let i = startPage; i <= endPage; i++) {{
                const isActive = i === currentSAGPage;
                paginationHTML += `<button onclick="changeSAGPage(${{i}})" style="padding: 0.5rem 1rem; background: ${{isActive ? '#495057' : '#f8f9fa'}}; color: ${{isActive ? 'white' : '#495057'}}; border: 1px solid #dee2e6; border-radius: 4px; cursor: pointer;">${{i}}</button>`;
            }}
            
            // Next button
            if (currentSAGPage < totalPages) {{
                paginationHTML += `<button onclick="changeSAGPage(${{currentSAGPage + 1}})" style="padding: 0.5rem 1rem; background: #667eea; color: white; border: none; border-radius: 4px; cursor: pointer;">Next</button>`;
            }}
            
            paginationHTML += `
                    </div>
                </div>
            `;
            
            paginationDiv.innerHTML = paginationHTML;
        }}

        // Change SAG page
        function changeSAGPage(page) {{
            currentSAGPage = page;
            displaySAGsTable();
        }}

        // Filter SAGs
        function filterSAGs() {{
            const searchTerm = document.getElementById('sag-search')?.value.toLowerCase() || '';
            const qualityFilter = document.getElementById('quality-filter')?.value || 'all';
            const taxonomyFilter = document.getElementById('taxonomy-filter')?.value || 'all';
            
            filteredSAGsData = allSAGsData.filter(item => {{
                const matchesSearch = !searchTerm || 
                    item.sag_id.toLowerCase().includes(searchTerm) ||
                    item.cluster_id.toLowerCase().includes(searchTerm) ||
                    item.taxonomy.toLowerCase().includes(searchTerm);
                
                const matchesQuality = qualityFilter === 'all' || 
                    item.quality.toLowerCase() === qualityFilter;
                
                let matchesTaxonomy = true;
                if (taxonomyFilter !== 'all') {{
                    if (taxonomyFilter === 'annotated') {{
                        matchesTaxonomy = item.taxonomy !== 'Unknown';
                    }} else if (taxonomyFilter === 'unknown') {{
                        matchesTaxonomy = item.taxonomy === 'Unknown';
                    }}
                }}
                
                return matchesSearch && matchesQuality && matchesTaxonomy;
            }});
            
            currentSAGPage = 1;
            displaySAGsTable();
        }}

        // Sort SAGs
        function sortSAGs() {{
            const sortBy = document.getElementById('sag-sort')?.value || 'sag_id';
            
            filteredSAGsData.sort((a, b) => {{
                switch (sortBy) {{
                    case 'completeness':
                        return b.completeness - a.completeness;
                    case 'contamination':
                        return a.contamination - b.contamination;
                    case 'quality':
                        const qualityOrder = {{ 'High': 3, 'Medium': 2, 'Low': 1 }};
                        return qualityOrder[b.quality] - qualityOrder[a.quality];
                    case 'cluster_id':
                        return a.cluster_id.localeCompare(b.cluster_id);
                    case 'sag_id':
                    default:
                        return a.sag_id.localeCompare(b.sag_id);
                }}
            }});
            
            currentSAGPage = 1;
            displaySAGsTable();
        }}

        // Download SAGs CSV
        function downloadSAGsCSV() {{
            const headers = ['SAG ID', 'Completeness (%)', 'Contamination (%)', 'Cluster ID', 'Quality', 'Taxonomy'];
            const csvContent = [
                headers.join(','),
                ...filteredSAGsData.map(item => [
                    item.sag_id,
                    item.completeness.toFixed(1),
                    item.contamination.toFixed(1),
                    item.cluster_id,
                    item.quality,
                    `"${{item.taxonomy}}"`
                ].join(','))
            ].join('\\n');
            
            const blob = new Blob([csvContent], {{ type: 'text/csv;charset=utf-8;' }});
            const link = document.createElement('a');
            const url = URL.createObjectURL(blob);
            link.setAttribute('href', url);
            link.setAttribute('download', 'single_sag_quality.csv');
            link.style.visibility = 'hidden';
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
        }}

        // Clustering Analysis functions
        let allClustersData = [];
        let filteredClustersData = [];
        let currentClusterPage = 1;
        const clusterItemsPerPage = 20;

        // Create clustering charts
        function createClusteringCharts() {{
            updateSingletonCount();
            createClusterSizeDistributionChart();
            createClusteringEfficiencyChart();
            createContaminationComparisonChart();
        }}

        // Create cluster size distribution chart with quality breakdown
        function createClusterSizeDistributionChart() {{
            const ctx = document.getElementById('cluster-size-distribution-chart');
            if (!ctx) return;

            // Collect data by size and quality
            const sizeQualityData = {{}};
            
            if (clusterData && clusterData.clusters) {{
                clusterData.clusters.forEach(cluster => {{
                    const size = cluster.cluster_size || 0;
                    const sizeKey = size > 10 ? '10+' : size.toString();
                    
                    // Determine cluster quality based on co-assembly results
                    const coResults = cluster.co_assembly_results || {{}};
                    const bestResult = coResults.optimized || coResults.round_1;
                    let quality = 'Low';
                    
                    if (bestResult) {{
                        const comp = bestResult.completeness || 0;
                        const cont = bestResult.contamination || 100;
                        
                        if (comp > 90 && cont < 5) {{
                            quality = 'High';
                        }} else if (comp > 50 && cont < 10) {{
                            quality = 'Medium';
                        }}
                    }}
                    
                    if (!sizeQualityData[sizeKey]) {{
                        sizeQualityData[sizeKey] = {{ High: 0, Medium: 0, Low: 0 }};
                    }}
                    sizeQualityData[sizeKey][quality]++;
                }});
            }}

            const sortedKeys = Object.keys(sizeQualityData).sort((a, b) => {{
                if (a === '10+') return 1;
                if (b === '10+') return -1;
                return parseInt(a) - parseInt(b);
            }});

            // Prepare data for stacked bar chart
            const highQualityData = sortedKeys.map(k => sizeQualityData[k]?.High || 0);
            const mediumQualityData = sortedKeys.map(k => sizeQualityData[k]?.Medium || 0);
            const lowQualityData = sortedKeys.map(k => sizeQualityData[k]?.Low || 0);

            new Chart(ctx, {{
                type: 'bar',
                data: {{
                    labels: sortedKeys.map(k => k + ' SAG' + (k === '1' ? '' : 's')),
                    datasets: [{{
                        label: 'High Quality Clusters',
                        data: highQualityData,
                        backgroundColor: '#28a745',
                        borderColor: '#1e7e34',
                        borderWidth: 1
                    }}, {{
                        label: 'Medium Quality Clusters',
                        data: mediumQualityData,
                        backgroundColor: '#ffc107',
                        borderColor: '#e0a800',
                        borderWidth: 1
                    }}, {{
                        label: 'Low Quality Clusters',
                        data: lowQualityData,
                        backgroundColor: '#dc3545',
                        borderColor: '#c82333',
                        borderWidth: 1
                    }}]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {{
                        legend: {{
                            display: true,
                            position: 'top',
                            labels: {{
                                padding: 20,
                                usePointStyle: true,
                                font: {{
                                    size: 12
                                }}
                            }}
                        }},
                        tooltip: {{
                            callbacks: {{
                                label: function(context) {{
                                    const datasetLabel = context.dataset.label;
                                    const value = context.parsed.y;
                                    return `${{datasetLabel}}: ${{value}}`;
                                }},
                                afterLabel: function(context) {{
                                    const dataIndex = context.dataIndex;
                                    const total = highQualityData[dataIndex] + mediumQualityData[dataIndex] + lowQualityData[dataIndex];
                                    const percentage = Math.round((context.parsed.y / total) * 100);
                                    return `(${{percentage}}% of size group)`;
                                }}
                            }}
                        }}
                    }},
                    scales: {{
                        x: {{
                            stacked: true,
                            title: {{
                                display: true,
                                text: 'Cluster Size'
                            }}
                        }},
                        y: {{
                            stacked: true,
                            beginAtZero: true,
                            title: {{
                                display: true,
                                text: 'Number of Clusters'
                            }}
                        }}
                    }}
                }}
            }});
        }}

        // Create clustering efficiency chart
        function createClusteringEfficiencyChart() {{
            const ctx = document.getElementById('clustering-efficiency-chart');
            if (!ctx) return;

            let totalSAGs = 0;
            let clusteredSAGs = 0;
            let singletonSAGs = 0;
            
            if (clusterData && clusterData.clusters) {{
                clusterData.clusters.forEach(cluster => {{
                    const size = cluster.cluster_size || 0;
                    totalSAGs += size;
                    if (size === 1) {{
                        singletonSAGs += size;
                    }} else {{
                        clusteredSAGs += size;
                    }}
                }});
            }}

            new Chart(ctx, {{
                type: 'doughnut',
                data: {{
                    labels: ['Successfully Clustered', 'Singletons'],
                    datasets: [{{
                        data: [clusteredSAGs, singletonSAGs],
                        backgroundColor: ['#28a745', '#ffc107'],
                        borderWidth: 2,
                        borderColor: '#fff'
                    }}]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {{
                        legend: {{
                            position: 'bottom',
                            labels: {{
                                padding: 20,
                                usePointStyle: true
                            }}
                        }},
                        tooltip: {{
                            callbacks: {{
                                label: function(context) {{
                                    const total = context.dataset.data.reduce((a, b) => a + b, 0);
                                    const percentage = Math.round((context.parsed / total) * 100);
                                    return `${{context.label}}: ${{context.parsed}} SAGs (${{percentage}}%)`;
                                }}
                            }}
                        }}
                    }}
                }}
            }});
        }}

        // Create contamination comparison chart
        function createContaminationComparisonChart() {{
            const ctx = document.getElementById('contamination-comparison-chart');
            if (!ctx) return;

            // Collect contamination data for clusters that have both round1 and optimized results
            const beforeData = [];
            const afterData = [];
            const labels = [];
            const reductions = [];
            
            if (clusterData && clusterData.clusters) {{
                clusterData.clusters.forEach(cluster => {{
                    const coResults = cluster.co_assembly_results || {{}};
                    const round1 = coResults.round_1;
                    const optimized = coResults.optimized;
                    
                    if (round1 && optimized && cluster.cluster_size > 1) {{
                        const beforeCont = round1.contamination || 0;
                        const afterCont = optimized.contamination || 0;
                        const reduction = beforeCont - afterCont; // Positive means contamination decreased (good)
                        
                        // Only show clusters with meaningful contamination changes (>0.2% reduction or any increase)
                        if (Math.abs(reduction) > 0.2) {{
                            beforeData.push(beforeCont);
                            afterData.push(afterCont);
                            labels.push(`Cluster ${{cluster.cluster_id}}`);
                            reductions.push(reduction);
                        }}
                    }}
                }});
            }}

            // If no meaningful changes, show a message
            if (labels.length === 0) {{
                const canvasCtx = ctx.getContext('2d');
                canvasCtx.font = '16px Arial';
                canvasCtx.fillStyle = '#6c757d';
                canvasCtx.textAlign = 'center';
                canvasCtx.fillText('No significant contamination changes detected', ctx.width/2, ctx.height/2);
                return;
            }}

            new Chart(ctx, {{
                type: 'bar',
                data: {{
                    labels: labels,
                    datasets: [{{
                        label: 'Before Optimization',
                        data: beforeData,
                        backgroundColor: '#dc3545',
                        borderColor: '#c82333',
                        borderWidth: 1
                    }}, {{
                        label: 'After Optimization',
                        data: afterData,
                        backgroundColor: '#28a745',
                        borderColor: '#1e7e34',
                        borderWidth: 1
                    }}]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {{
                        legend: {{
                            position: 'top',
                            labels: {{
                                padding: 20,
                                usePointStyle: true
                            }}
                        }},
                        tooltip: {{
                            callbacks: {{
                                afterLabel: function(context) {{
                                    const index = context.dataIndex;
                                    const reduction = reductions[index];
                                    if (reduction > 0) {{
                                        return `Contamination reduced by ${{reduction.toFixed(1)}}%`;
                                    }} else {{
                                        return `Contamination increased by ${{Math.abs(reduction).toFixed(1)}}%`;
                                    }}
                                }}
                            }}
                        }}
                    }},
                    scales: {{
                        x: {{
                            title: {{
                                display: true,
                                text: 'Clusters with Contamination Changes'
                            }}
                        }},
                        y: {{
                            title: {{
                                display: true,
                                text: 'Contamination (%)'
                            }},
                            beginAtZero: true,
                            max: function(context) {{
                                const maxBefore = Math.max(...beforeData);
                                const maxAfter = Math.max(...afterData);
                                return Math.max(maxBefore, maxAfter) * 1.1;
                            }}
                        }}
                    }}
                }}
            }});
        }}

        // Update singleton clusters count
        function updateSingletonCount() {{
            let singletonCount = 0;
            if (clusterData && clusterData.clusters) {{
                clusterData.clusters.forEach(cluster => {{
                    if ((cluster.cluster_size || 0) === 1) {{
                        singletonCount++;
                    }}
                }});
            }}
            const singletonElement = document.getElementById('singleton-clusters');
            if (singletonElement) {{
                singletonElement.textContent = singletonCount;
            }}
        }}

        // Populate clusters table
        function populateClustersTable() {{
            allClustersData = [];
            
            if (clusterData && clusterData.clusters) {{
                clusterData.clusters.forEach(cluster => {{
                    const coResults = cluster.co_assembly_results || {{}};
                    const round1 = coResults.round_1;
                    const optimized = coResults.optimized;
                    const bestResult = optimized || round1;
                    
                    // Get taxonomy - use full classification string
                    let taxonomy = 'Unknown';
                    if (bestResult && bestResult.gtdbtk_classification && bestResult.gtdbtk_classification.classification_string) {{
                        taxonomy = bestResult.gtdbtk_classification.classification_string;
                    }}
                    
                    // Calculate detailed improvement
                    let improvement = 'No optimization';
                    if (round1 && optimized) {{
                        const round1Comp = round1.completeness || 0;
                        const round1Cont = round1.contamination || 0;
                        const optComp = optimized.completeness || 0;
                        const optCont = optimized.contamination || 0;
                        
                        const compChange = optComp - round1Comp;
                        const contChange = round1Cont - optCont; // Positive means contamination decreased (good)
                        
                        let improvements = [];
                        if (Math.abs(compChange) >= 0.1) {{
                            improvements.push(`Completeness ${{compChange > 0 ? '+' : ''}}${{compChange.toFixed(1)}}%`);
                        }}
                        if (Math.abs(contChange) >= 0.1) {{
                            improvements.push(`Contamination ${{contChange > 0 ? '-' : '+'}}${{Math.abs(contChange).toFixed(1)}}%`);
                        }}
                        
                        if (improvements.length > 0) {{
                            improvement = improvements.join('; ');
                        }} else {{
                            improvement = 'Minimal change';
                        }}
                    }}
                    
                    // Get member list
                    const members = cluster.members || [];
                    const memberList = members.map(m => m.sag_id || 'Unknown').join(', ');
                    
                    allClustersData.push({{
                        cluster_id: cluster.cluster_id || 'Unknown',
                        size: cluster.cluster_size || 0,
                        members: memberList,
                        completeness: bestResult ? (bestResult.completeness || 0) : 0,
                        contamination: bestResult ? (bestResult.contamination || 0) : 0,
                        improvement: improvement,
                        taxonomy: taxonomy
                    }});
                }});
            }}
            
            filteredClustersData = [...allClustersData];
            displayClustersTable();
        }}

        // Display clusters table with pagination
        function displayClustersTable() {{
            const tbody = document.getElementById('cluster-table-body');
            if (!tbody) return;
            
            const startIndex = (currentClusterPage - 1) * clusterItemsPerPage;
            const endIndex = startIndex + clusterItemsPerPage;
            const pageData = filteredClustersData.slice(startIndex, endIndex);
            
            tbody.innerHTML = pageData.map(item => `
                <tr>
                    <td>${{item.cluster_id}}</td>
                    <td>${{item.size}}</td>
                    <td title="${{item.members}}">${{item.members.length > 50 ? item.members.substring(0, 50) + '...' : item.members}}</td>
                    <td>${{item.completeness.toFixed(1)}}</td>
                    <td>${{item.contamination.toFixed(1)}}</td>
                    <td style="font-size: 0.9rem;">${{item.improvement}}</td>
                    <td style="font-style: italic; max-width: 300px; word-wrap: break-word; font-size: 0.85rem;" title="${{item.taxonomy}}">${{item.taxonomy}}</td>
                </tr>
            `).join('');
            
            updateClusterPagination();
        }}

        // Update cluster pagination controls
        function updateClusterPagination() {{
            const paginationDiv = document.getElementById('cluster-pagination');
            if (!paginationDiv) return;
            
            const totalPages = Math.ceil(filteredClustersData.length / clusterItemsPerPage);
            const totalItems = filteredClustersData.length;
            const startItem = (currentClusterPage - 1) * clusterItemsPerPage + 1;
            const endItem = Math.min(currentClusterPage * clusterItemsPerPage, totalItems);
            
            let paginationHTML = `
                <div style="display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 1rem;">
                    <div>Showing ${{startItem}}-${{endItem}} of ${{totalItems}} clusters</div>
                    <div style="display: flex; gap: 0.5rem; align-items: center;">
            `;
            
            // Previous button
            if (currentClusterPage > 1) {{
                paginationHTML += `<button onclick="changeClusterPage(${{currentClusterPage - 1}})" style="padding: 0.5rem 1rem; background: #667eea; color: white; border: none; border-radius: 4px; cursor: pointer;">Previous</button>`;
            }}
            
            // Page numbers
            const maxVisiblePages = 5;
            let startPage = Math.max(1, currentClusterPage - Math.floor(maxVisiblePages / 2));
            let endPage = Math.min(totalPages, startPage + maxVisiblePages - 1);
            
            if (endPage - startPage + 1 < maxVisiblePages) {{
                startPage = Math.max(1, endPage - maxVisiblePages + 1);
            }}
            
            for (let i = startPage; i <= endPage; i++) {{
                const isActive = i === currentClusterPage;
                paginationHTML += `<button onclick="changeClusterPage(${{i}})" style="padding: 0.5rem 1rem; background: ${{isActive ? '#495057' : '#f8f9fa'}}; color: ${{isActive ? 'white' : '#495057'}}; border: 1px solid #dee2e6; border-radius: 4px; cursor: pointer;">${{i}}</button>`;
            }}
            
            // Next button
            if (currentClusterPage < totalPages) {{
                paginationHTML += `<button onclick="changeClusterPage(${{currentClusterPage + 1}})" style="padding: 0.5rem 1rem; background: #667eea; color: white; border: none; border-radius: 4px; cursor: pointer;">Next</button>`;
            }}
            
            paginationHTML += `
                    </div>
                </div>
            `;
            
            paginationDiv.innerHTML = paginationHTML;
        }}

        // Change cluster page
        function changeClusterPage(page) {{
            currentClusterPage = page;
            displayClustersTable();
        }}

        // Filter clusters
        function filterClusters() {{
            const searchTerm = document.getElementById('cluster-search')?.value.toLowerCase() || '';
            const sizeFilter = document.getElementById('cluster-size-filter')?.value || 'all';
            const taxonomyFilter = document.getElementById('cluster-taxonomy-filter')?.value || 'all';
            
            filteredClustersData = allClustersData.filter(item => {{
                const matchesSearch = !searchTerm || 
                    item.cluster_id.toLowerCase().includes(searchTerm) ||
                    item.members.toLowerCase().includes(searchTerm) ||
                    item.taxonomy.toLowerCase().includes(searchTerm) ||
                    item.improvement.toLowerCase().includes(searchTerm);
                
                let matchesSize = true;
                if (sizeFilter !== 'all') {{
                    if (sizeFilter === '1') {{
                        matchesSize = item.size === 1;
                    }} else if (sizeFilter === '2') {{
                        matchesSize = item.size === 2;
                    }} else if (sizeFilter === '3-5') {{
                        matchesSize = item.size >= 3 && item.size <= 5;
                    }} else if (sizeFilter === '6+') {{
                        matchesSize = item.size >= 6;
                    }}
                }}
                
                let matchesTaxonomy = true;
                if (taxonomyFilter !== 'all') {{
                    if (taxonomyFilter === 'unknown') {{
                        matchesTaxonomy = item.taxonomy === 'Unknown' || item.taxonomy === 'N/A' || item.taxonomy === '';
                    }} else if (taxonomyFilter === 'annotated') {{
                        matchesTaxonomy = item.taxonomy !== 'Unknown' && item.taxonomy !== 'N/A' && item.taxonomy !== '';
                    }}
                }}
                
                return matchesSearch && matchesSize && matchesTaxonomy;
            }});
            
            currentClusterPage = 1;
            displayClustersTable();
        }}

        // Sort clusters
        function sortClusters() {{
            const sortBy = document.getElementById('cluster-sort')?.value || 'cluster_id';
            
            filteredClustersData.sort((a, b) => {{
                switch (sortBy) {{
                    case 'size':
                        return b.size - a.size;
                    case 'completeness':
                        return b.completeness - a.completeness;
                    case 'contamination':
                        return a.contamination - b.contamination;
                    case 'improvement':
                        // Parse improvement values for sorting
                        const aImpr = a.improvement === 'N/A' ? -999 : parseFloat(a.improvement.replace('%', ''));
                        const bImpr = b.improvement === 'N/A' ? -999 : parseFloat(b.improvement.replace('%', ''));
                        return bImpr - aImpr;
                    case 'cluster_id':
                    default:
                        return a.cluster_id.localeCompare(b.cluster_id);
                }}
            }});
            
            currentClusterPage = 1;
            displayClustersTable();
        }}

        // Download clusters CSV
        function downloadClustersCSV() {{
            const headers = ['Cluster ID', 'Size', 'Members', 'Completeness (%)', 'Contamination (%)', 'Quality Improvement', 'Taxonomy'];
            const csvContent = [
                headers.join(','),
                ...filteredClustersData.map(item => [
                    item.cluster_id,
                    item.size,
                    `"${{item.members}}"`,
                    item.completeness.toFixed(1),
                    item.contamination.toFixed(1),
                    `"${{item.improvement}}"`,
                    `"${{item.taxonomy}}"`
                ].join(','))
            ].join('\\n');
            
            const blob = new Blob([csvContent], {{ type: 'text/csv;charset=utf-8;' }});
            const link = document.createElement('a');
            const url = URL.createObjectURL(blob);
            link.setAttribute('href', url);
            link.setAttribute('download', 'clustering_analysis.csv');
            link.style.visibility = 'hidden';
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
        }}

        // Initialize the page
        document.addEventListener('DOMContentLoaded', function() {{
            createQualityChart();
        }});
    </script>
</body>
</html>"""

    # Write HTML file
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    return True


def main():
    parser = argparse.ArgumentParser(description='Create HTML report with embedded data')
    parser.add_argument('json_file', help='Input cluster data JSON file')
    parser.add_argument('-o', '--output', default='embedded_report.html', help='Output HTML file')
    parser.add_argument('-t', '--title', default='SAG Analysis Report', help='Report title')
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.json_file):
        print(f"Error: Input file '{args.json_file}' not found!")
        sys.exit(1)
    
    # Load data
    print(f"Loading data from {args.json_file}...")
    data = load_cluster_data(args.json_file)
    if not data:
        sys.exit(1)
    
    # Generate report
    print(f"Creating embedded HTML report...")
    success = create_embedded_html_report(data, args.title, args.output)
    
    if success:
        print(f"HTML report created successfully: {args.output}")
        
        # Show summary
        stats = calculate_summary_stats(data)
        if stats:
            print(f"\nReport Summary:")
            print(f"  Total SAGs: {stats['total_sags']}")
            print(f"  High and Medium Quality CoSAGs: {stats['high_medium_quality_cosags']}")
            print(f"  Total Clusters: {stats['total_clusters']}")
            print(f"  Optimized Clusters: {stats['optimized_clusters']}")
            print(f"  Final Quality Genomes: {stats['final_quality_genomes']}")
        
        print(f"\nYou can now open {args.output} directly in your browser!")
    else:
        print("Failed to create report!")
        sys.exit(1)


if __name__ == "__main__":
    main()