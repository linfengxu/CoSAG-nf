# linfengxu/cosag-nf: Usage

## Introduction

**cosag-nf** processes paired-end single amplified genomes (SAGs): individual SPAdes assembly, Sourmash clustering, co-assembly, optional TNF optimization, GTDB-Tk classification, and an interactive HTML report.

For installation, containers, and databases, see the [main README](../README.md).

## Input data prerequisites

CoSAG-nf expects **preprocessed paired-end FASTQ** files (validated on Illumina short reads; long-read and single-end data are not supported). Raw reads should be quality-filtered and, when applicable, host-decontaminated **before** building the samplesheet and launching the pipeline.

| Step | Purpose | Recommended tools |
| ---- | ------- | ----------------- |
| Quality control | Adapter trimming, length/quality filtering | [fastp](https://github.com/OpenGene/fastp), [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) |
| Host removal | Remove human or other host reads (e.g. oral microbiome samples) | [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/) (manual workflow), [KneadData](https://huttenhower.sph.harvard.edu/kneaddata/) |

Point `forwardReads` and `reverseReads` in the samplesheet to the **final cleaned FASTQs**. CoSAG-nf does not run read QC or host decontamination internally.

## Pipeline parameters

Full parameter documentation is maintained in [`nextflow_schema.json`](../nextflow_schema.json), grouped by topic:

| Schema section                   | Examples                                                                              |
| -------------------------------- | ------------------------------------------------------------------------------------- |
| Input/output options             | `--input`, `--outdir`, `keep_intermediate`                                            |
| Assembly options                 | `--spades_kmers`, `--spades_sc`, `--spades_careful`                                   |
| MinHash clustering options       | `--sourmash_ksize`, `--cluster_threshold`, `--cluster_min_size`                       |
| Quality filtering options        | `--checkm2_db`, `--cosag_rounds`, `--round`, `--cosag_optimize`, `--min_completeness` |
| TNF optimization options         | `--tnf_max_iterations`, `--tnf_min_sags`                                              |
| Taxonomic classification options | `--gtdb_database`                                                                     |
| Output structure options         | `--out_individual_assemblies`, … `--out_logs_diagnostics`                             |
| Max resource options             | `--max_cpus`, `--max_memory`, `--max_time`                                            |

Defaults live in [`nextflow.config`](../nextflow.config). Copy [`assets/params.example.yaml`](../assets/params.example.yaml) as a starting point for `-params-file`.

Browse parameters interactively: [nf-core launch](https://nf-co.re/launch) (select `linfengxu/cosag-nf`) or open `nextflow_schema.json` in the repository.

> [!NOTE]
> At startup the pipeline prints a parameter summary and validates inputs when `--validate_params` is true (default).

### Commonly confused parameters

| Parameter                                                        | Meaning                                                                              |
| ---------------------------------------------------------------- | ------------------------------------------------------------------------------------ |
| `--round 2`                                                      | Run Round 1, then a **separate** Round 2 workflow; outputs under `<outdir>/round2/`. |
| `--min_completeness` / `--max_contamination`                     | Thresholds for selecting final CoSAG contigs.                                        |
| `--cosag_opt_min_completeness` / `--cosag_opt_max_contamination` | When to **trigger** TNF optimization on a cluster.                                   |
| `--target_completeness` / `--target_contamination`               | Targets used **inside** the TNF optimizer iterations.                                |

### Hierarchical clustering linkage

SAG grouping uses SciPy hierarchical clustering on the Sourmash distance matrix. Supported linkage methods (`--cluster_linkage_method`): `complete` (default), `average`, `single`, `ward`. This is a **CLI parameter only** — no extra linkage input file.

| Method     | Summary                                                                                   |
| ---------- | ----------------------------------------------------------------------------------------- |
| `complete` | Maximizes minimum intra-cluster similarity; conservative default for co-assembly clusters |
| `average`  | Average inter-cluster distance                                                            |
| `single`   | Minimum distance; watch for chaining into oversized clusters                              |
| `ward`     | Variance-based merges; validate on your data if used here                                 |

Cut the tree with `--cluster_criterion` (default `inconsistent`) and `--cluster_threshold` (default `0.95`). Re-tune the threshold when changing linkage. See [`nextflow_schema.json`](../nextflow_schema.json) (`minhash_clustering_options`) for full `help_text`, and the [README linkage section](../README.md#hierarchical-clustering-linkage-methods).

## Samplesheet input

Create a tab- or comma-separated samplesheet and pass it with `--input`. The file is validated against [`assets/schema_input.json`](../assets/schema_input.json).

```bash
--input /path/to/samples.tsv
```

### Required columns

| Column         | Description                                                                                                            |
| -------------- | ---------------------------------------------------------------------------------------------------------------------- |
| `sampleID`     | Unique SAG identifier. Use the same ID on multiple rows to merge sequencing lanes (reads are concatenated). No spaces. |
| `forwardReads` | Absolute path to read 1 (`.fastq`, `.fq`, `.fastq.gz`, or `.fq.gz`).                                                   |
| `reverseReads` | Absolute path to read 2 (same extensions). Required by the input schema for all rows.                                  |

Example (TSV):

```tsv
sampleID	forwardReads	reverseReads
SAG001	/path/to/SAG001_R1.fastq.gz	/path/to/SAG001_R2.fastq.gz
SAG002	/path/to/SAG002_R1.fastq.gz	/path/to/SAG002_R2.fastq.gz
```

### Multiple lanes per sample

Repeat the same `sampleID` on multiple rows; the pipeline merges FASTQs before assembly:

```tsv
sampleID	forwardReads	reverseReads
SAG001	/path/to/SAG001_L002_R1.fastq.gz	/path/to/SAG001_L002_R2.fastq.gz
SAG001	/path/to/SAG001_L003_R1.fastq.gz	/path/to/SAG001_L003_R2.fastq.gz
```

An example file is provided: [`assets/samplesheet.csv`](../assets/samplesheet.csv).

## Running the pipeline

### Minimal command

```bash
nextflow run linfengxu/cosag-nf \
    -profile singularity \
    --input samples.tsv \
    --outdir results \
    --checkm2_db /path/to/checkm2_database/uniref100.KO.1.dmnd \
    --gtdb_database /path/to/gtdbtk_r220_data
```

From a local clone:

```bash
nextflow run main.nf -profile singularity --input samples.tsv --outdir results ...
```

### Multi-round and Round 2

**Internal co-assembly rounds** (`--cosag_rounds` 2–3): extra co-assembly iterations inside Round 1; JSON under `04_co_assemblies/cluster_json/round2/` etc.

**Two-round sourmash clustering** (`--round 2`): after Round 1 completes, a second MinHash pass clusters **Round 1 co-assembly contigs** (cluster representatives), then co-assembly and reporting run again under `<outdir>/round2/`. Useful when uneven SAG coverage weakens first-pass sketch overlap. Distances use **1 − Jaccard similarity**. For SAG-level clustering, try `--sourmash_ksize 31` (default `51` was strict in an oral microbiome benchmark).

```bash
# Internal co-assembly rounds 2–3 (same outdir, round-specific cluster JSON)
nextflow run linfengxu/CoSAG-nf -profile singularity \
    --input samples.tsv --outdir results --cosag_rounds 2

# Round 1 + second sourmash pass on co-assembly contigs (results/round2/)
nextflow run linfengxu/CoSAG-nf -profile singularity \
    --input samples.tsv --outdir results --round 2 --sourmash_ksize 31
```

### Parameter file

```bash
nextflow run linfengxu/cosag-nf -profile singularity -params-file assets/params.example.yaml
```

> [!WARNING]
> Do not pass pipeline parameters with `-c`. Use CLI flags or `-params-file` only. Use `-c` for resource tuning or institutional config snippets.

### Working directory

```text
work/           # Nextflow work directory
<OUTDIR>/       # Published results (--outdir)
.nextflow_log   # Nextflow log
```

See [Output](output.md) for the structure under `<OUTDIR>/`.

### Updating and versioning

```bash
nextflow pull linfengxu/cosag-nf
nextflow run linfengxu/cosag-nf -r <version> ...
```

## Core Nextflow arguments

> [!NOTE]
> Nextflow options use a **single** hyphen; pipeline parameters use a **double** hyphen.

### `-profile`

Recommended: `singularity` or `docker` (see [containers/README.md](../containers/README.md)).

| Profile       | Description                                  |
| ------------- | -------------------------------------------- |
| `test`        | Small test dataset; minimal extra parameters |
| `docker`      | Run with Docker                              |
| `singularity` | Run with Singularity / Apptainer             |
| `conda`       | Conda environments (fallback only)           |

Multiple profiles can be combined, e.g. `-profile test,singularity` (order matters).

Institutional configs are loaded from [nf-core/configs](https://github.com/nf-core/configs) when available.

### `-resume`

Resume a previous run; reuse cached tasks with identical inputs.

```bash
nextflow run linfengxu/cosag-nf -profile singularity -resume
```

### `-c`

Path to an extra Nextflow config file for resources or cluster settings — **not** for overriding `params.*`.

## Custom configuration

### Resource limits

Adjust global caps in the schema / config: `--max_cpus`, `--max_memory`, `--max_time`. Per-process resources are in [`conf/base.config`](../conf/base.config) and can be overridden with `-c` (see [nf-core configuration docs](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources)).

### Containers

CoSAG-specific images are documented in [`containers/README.md`](../containers/README.md). Other tools may still use local `.sif` paths in `modules/local/*/main.nf`.

## Running in the background

```bash
nextflow run linfengxu/cosag-nf -profile singularity -bg ...
```

Or use `screen` / `tmux`, or submit Nextflow itself as a cluster job.

## Nextflow JVM memory

```bash
export NXF_OPTS='-Xms1g -Xmx4g'
```

Add to `~/.bashrc` if Nextflow uses excessive memory on the login node.
