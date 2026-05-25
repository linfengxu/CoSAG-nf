# linfengxu/cosag-nf: Output

## Introduction

This document describes the output produced by the **cosag-nf** pipeline. The pipeline assembles single amplified genomes (SAGs), clusters them by MinHash similarity, performs iterative co-assembly, applies TNF-based SAG subset optimization, and produces taxonomically annotated CoSAG (co-assembled SAG) genomes with an interactive HTML summary report.

All result paths below are relative to the top-level output directory specified with `--outdir`. Directory names `01_individual_assemblies` through `07_logs_and_diagnostics` can be customized via the corresponding `out_*` parameters in `nextflow.config`.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data in the following stages:

1. [Individual SAG assembly](#individual-sag-assembly) — SPAdes assembly, CheckM2 QC, contamination filtering, GTDB-Tk classification
2. [Similarity analysis](#similarity-analysis) — Sourmash MinHash sketching and pairwise similarity matrix
3. [Clustering analysis](#clustering-analysis) — Hierarchical clustering and cluster JSON generation
4. [Co-assemblies](#co-assemblies) — Merged-read co-assembly, CheckM2, optional TNF optimization, cluster JSON updates
5. [Taxonomic classification](#taxonomic-classification) — GTDB-Tk on individual SAGs and final CoSAGs
6. [Final results](#final-results) — Integrated cluster JSON and HTML report
7. [Logs and diagnostics](#logs-and-diagnostics) — Filtering logs, CheckM2 run logs, parsed cluster tables
8. [Pipeline information](#pipeline-information) — Nextflow execution reports and run metadata

When `--round 2` is set, the pipeline first completes Round 1 (steps 1–8 above), then runs an additional **Round 2** workflow. Round 2 outputs are written under `<OUTDIR>/round2/` with the same numbered subdirectory layout (`02_similarity_analysis` … `06_final_results`).

When `--cosag_rounds` is greater than 1, additional internal co-assembly rounds (Round 2/3) are executed within the main workflow. Per-round cluster JSON files are stored under `04_co_assemblies/cluster_json/round{N}/`.

## Output directory layout

A typical successful run produces the following structure:

```
<OUTDIR>/
├── 01_individual_assemblies/
├── 02_similarity_analysis/
├── 03_clustering_analysis/
├── 04_co_assemblies/
├── 05_taxonomic_classification/
├── 06_final_results/
├── 07_logs_and_diagnostics/
├── pipeline_info/
└── round2/                          # only when --round 2
    ├── 02_similarity_analysis/
    ├── 03_clustering_analysis/
    ├── 04_co_assemblies/
    ├── 05_taxonomic_classification/
    ├── 06_final_results/
    └── 07_logs_and_diagnostics/
```

### Individual SAG assembly

<details markdown="1">
<summary>Output files</summary>

- `01_individual_assemblies/per_sample/<SAMPLE_ID>/`
  - `contigs/<SAMPLE_ID>_contigs.fasta`: SPAdes assembled contigs for each SAG.
  - `logs/<SAMPLE_ID>_spades.log`: SPAdes assembly log.
  - `graphs/<SAMPLE_ID>_assembly_graph.gfa`: Assembly graph (published when available).
- `01_individual_assemblies/checkm2/checkm2_results/quality_report.tsv`: Batch CheckM2 quality metrics for all individual SAGs.

</details>

Individual SAGs are assembled with SPAdes in single-cell mode, evaluated with [CheckM2](https://github.com/CheckM-CheckM2/CheckM2), and filtered by contamination (`FILTER_SAGS`, default threshold: contamination > 10%). SAGs passing the filter proceed to downstream clustering and co-assembly.

### Similarity analysis

<details markdown="1">
<summary>Output files</summary>

- `02_similarity_analysis/signatures/`
  - `<SAMPLE_ID>_k51.sig`: Sourmash MinHash signatures (k-mer size controlled by `--sourmash_ksize`).
- `02_similarity_analysis/matrix/`
  - `sourmash_similarity_matrix.tsv`: Pairwise similarity matrix between SAG sketches.
  - `sourmash_similarity_report.txt`: Summary statistics of the similarity matrix.
  - `sourmash_comparison_details.csv`: Per-comparison details from Sourmash.
  - `sourmash_compare.log`: Sourmash compare process log.
  - `failed_comparisons.txt`: Comparisons that failed (if any).
- `02_similarity_analysis/quality_filter/` *(optional, when `--sourmash_enable_quality_filter` is true)*
  - Quality-filtered sketch lists and statistics.

</details>

[Sourmash](https://sourmash.readthedocs.io/) compares MinHash sketches to quantify genomic similarity between SAGs, forming the basis for hierarchical clustering.

### Clustering analysis

<details markdown="1">
<summary>Output files</summary>

- `03_clustering_analysis/matrix/`
  - `similarity_stats.txt`: Distance/similarity statistics used for clustering.
- `03_clustering_analysis/results/`
  - Cluster assignment tables (`*.tsv`, `*.txt`).
  - `dendrogram.png`: Hierarchical clustering dendrogram (when `--cluster_generate_dendrogram` is true).
- `03_clustering_analysis/summary/`
  - `final_clustering_summary.txt`: Overview of clustering results.
  - `cluster_statistics.tsv`: Per-cluster size and quality statistics.
  - `high_quality_clusters.tsv`: Clusters meeting quality criteria for co-assembly.
- `03_clustering_analysis/cluster_json/round1/`
  - `round1_clusters.json`: Cluster definitions with member SAGs and merged read paths for Round 1 co-assembly.

</details>

Clusters are built with hierarchical clustering (linkage method and threshold controlled by `--cluster_linkage_method` and `--cluster_threshold`). Each cluster in `round1_clusters.json` lists member SAG IDs and paths to merged paired-end reads used for co-assembly.

### Co-assemblies

<details markdown="1">
<summary>Output files</summary>

- `04_co_assemblies/coassembly/<CLUSTER_ID>/`
  - `contigs/<CLUSTER_ID>_contigs.fasta`: SPAdes co-assembly contigs per cluster (e.g. `R1_C602_contigs.fasta`).
  - `logs/<CLUSTER_ID>_spades.log`: Co-assembly SPAdes log.
- `04_co_assemblies/checkm2/checkm2_results/quality_report.tsv`: CheckM2 metrics for all co-assembled genomes.
- `04_co_assemblies/cluster_json/`
  - `round1/updated_clusters.json`: Cluster JSON after Round 1 co-assembly (includes `coassembly_contigs` and `checkm2_coassembly`).
  - `round2/updated_clusters.json`, `round3/updated_clusters.json`: Per-round updated JSON when `--cosag_rounds` ≥ 2 or 3.
  - `latest/updated_clusters.json`: Most recent cluster JSON snapshot.
- `04_co_assemblies/tnf_optimization/results/<CLUSTER_ID>/`
  - `*_tnf_optimized.json`: TNF optimization results for clusters meeting optimization criteria (completeness > `--cosag_opt_min_completeness`, contamination > `--cosag_opt_max_contamination`, with `--cosag_optimize` enabled).
- `04_co_assemblies/tnf_optimization/logs/<CLUSTER_ID>/`
  - `*_sag_opt.log`: TNF optimizer logs.
- `04_co_assemblies/selected_cosags/`
  - `fasta/selected_fastas/cluster_<CLUSTER_ID>.fasta`: Final CoSAG contigs passing completeness/contamination thresholds.
  - `metadata/selected_fastas/selected_contigs.tsv`: Mapping of selected contigs to clusters.
  - `barrnap/`: rRNA annotation (GFF and logs) from barrnap when `--run_cosag_rrna_annotation` is true.
- `04_co_assemblies/filtered_clusters/json/`: JSON records for clusters removed due to high contamination (when filtering is applied).
- `04_co_assemblies/merged_reads/<CLUSTER_ID>/`: Logs from read merging prior to co-assembly.

</details>

Co-assembly merges reads from all SAGs within a cluster and assembles them with SPAdes. For clusters with high completeness but elevated contamination, the TNF optimizer selects a subset of member SAGs to improve genome quality without re-assembly. Cluster IDs follow the pattern `R{N}_C<id>` (e.g. `R1_C602` for Round 1, `R2_C730` for Round 2).

### Taxonomic classification

<details markdown="1">
<summary>Output files</summary>

- `05_taxonomic_classification/individual_sags/gtdb_results/gtdb_results/`
  - Standard [GTDB-Tk](https://ecogenomics.github.io/GTDBTk/) output for passing individual SAGs, including `gtdbtk.bac120.summary.tsv` and `gtdbtk.ar53.summary.tsv`.
- `05_taxonomic_classification/co_assemblies/gtdb_results/gtdb_results/`
  - GTDB-Tk classification results for selected CoSAG genomes (same file layout as above).

</details>

Taxonomic assignments from GTDB-Tk are merged into the final cluster JSON (see below). Classification uses the database specified by `--gtdb_database`.

### Final results

<details markdown="1">
<summary>Output files</summary>

- `06_final_results/cluster_json/`
  - `cluster_data_tnf.json`: Cluster metadata with TNF optimization results integrated.
  - `cluster_data_gtdbtk.json`: **Primary deliverable** — full cluster JSON with CheckM2 metrics, GTDB-Tk classifications (SAG members and CoSAGs), barrnap rRNA counts, and HQ MAG assessment flags (`hq_mag_final` when `--hq_mag_final_assessment` is true).
- `06_final_results/report/`
  - `cosag_report.html`: Standalone interactive HTML report with embedded cluster data. Open in a web browser to explore per-cluster assembly quality, taxonomy, and optimization history.

</details>

The `cluster_data_gtdbtk.json` file is the central structured output. Each cluster entry includes member SAGs, co-assembly results per round, CheckM2 metrics, GTDB-Tk classification strings, and optional HQ MAG status. Downstream analyses should start from this file.

### Logs and diagnostics

<details markdown="1">
<summary>Output files</summary>

- `07_logs_and_diagnostics/filter_sags/`
  - `sags_fail.tsv`: SAGs failing contamination filter.
  - `filter_stats.json`: Summary counts of passed/failed SAGs.
- `07_logs_and_diagnostics/checkm2/`
  - `individual_sags/checkm2_results.log`: CheckM2 run log for individual SAGs.
  - `coassemblies/checkm2_results.log`: CheckM2 run log for co-assemblies.
- `07_logs_and_diagnostics/filtering/high_contamination/`
  - `filter_contam.log`: Log from high-contamination cluster filtering.
- `07_logs_and_diagnostics/cluster_json/parse/`
  - `clusters.tsv`: Flattened cluster membership table parsed from cluster JSON.

</details>

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - `execution_report.html`, `execution_timeline.html`, `execution_trace.txt`: Nextflow execution reports.
  - `pipeline_dag.html`: Workflow DAG.
  - `params_<timestamp>.json`: Parameters used for the run.
  - `software_versions.yml`: Software versions (when version reporting is enabled).
  - `pipeline_report.html`, `pipeline_report.txt`: Email summary reports (only when `--email` or `--email_on_fail` is set).

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) tracing reports help troubleshoot failed steps and review resource usage.

## Round 2 outputs

When the pipeline is run with `--round 2`, Round 1 completes first and Round 2 re-clusters Round 1 co-assembly contigs (excluding clusters with contamination above `--max_contamination`). Round 2 results are published under `<OUTDIR>/round2/`:

| Round 1 path | Round 2 equivalent |
|---|---|
| `03_clustering_analysis/cluster_json/round1/` | `round2/03_clustering_analysis/cluster_json/round2/round2_clusters.json` |
| `04_co_assemblies/cluster_json/round1/` | `round2/04_co_assemblies/cluster_json/round2/updated_clusters.json` |
| `04_co_assemblies/coassembly/R1_*` | `round2/04_co_assemblies/coassembly/R2_*` |
| `06_final_results/cluster_json/cluster_data_gtdbtk.json` | `round2/06_final_results/cluster_json/cluster_data_gtdbtk.json` |
| `06_final_results/report/cosag_report.html` | `round2/06_final_results/report/cosag_report.html` |

Round 1 outputs remain in the top-level `<OUTDIR>/` directories and are not overwritten by Round 2.

## Key parameters affecting output

| Parameter | Default | Effect on output |
|---|---|---|
| `--cosag_rounds` | `1` | Number of internal co-assembly rounds (1–3); adds `round2/` and `round3/` under `04_co_assemblies/cluster_json/` |
| `--round` | `1` | Set to `2` to enable the standalone Round 2 workflow and `round2/` output subdirectory |
| `--cosag_optimize` | `true` | Enables TNF optimization outputs under `tnf_optimization/` |
| `--min_completeness` / `--max_contamination` | `50` / `10` | Thresholds for selecting CoSAG contigs in `selected_cosags/` |
| `--hq_mag_final_assessment` | `true` | Adds HQ MAG flags to `cluster_data_gtdbtk.json` |
| `--keep_intermediate` | `false` | Retains additional intermediate files in the work directory (not published by default) |
| `--publish_dir_mode` | `copy` | How files are written to the output directory (`copy`, `symlink`, or `rellink`) |

Output subdirectory names can be overridden individually, for example `--out_co_assemblies my_coassembly_dir`.
