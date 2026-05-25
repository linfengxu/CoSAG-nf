# linfengxu/CoSAG-nf: Documentation

The [linfengxu/CoSAG-nf](https://github.com/linfengxu/CoSAG-nf) documentation is split into the following pages:

- [Usage](usage.md)
  - How to run the pipeline, samplesheet format, and where to find parameter help.
- [Output](output.md)
  - Output directory layout, key files, and how to interpret results.

## Pipeline parameters

Command-line parameters (`--input`, `--outdir`, `--cosag_rounds`, etc.) are **not** listed in full in these markdown files. They are defined and documented in:

| Resource                                                      | Description                                                                                                                                                                    |
| ------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| [`nextflow_schema.json`](../nextflow_schema.json)             | Authoritative parameter definitions: `description`, `help_text`, defaults, types, and validation rules. Used by [nf-core launch](https://nf-co.re/launch) and Seqera Platform. |
| [`nextflow.config`](../nextflow.config)                       | Default values for each `params.*` entry (must stay in sync with the schema).                                                                                                  |
| [`assets/params.example.yaml`](../assets/params.example.yaml) | Example parameter file for `-params-file`.                                                                                                                                     |

At run time the pipeline validates parameters against the schema (`UTILS_NFSCHEMA_PLUGIN`) and prints a summary of values that differ from defaults.

To edit parameter help text, update `description` / `help_text` in `nextflow_schema.json` under the relevant `$defs` section (for example `quality_filtering_options`, `minhash_clustering_options`).

```bash
# Optional: validate schema against nextflow.config (requires nf-core tools)
nf-core pipelines schema validate
```

## Input samplesheet schema

The format of `--input` is validated separately via [`assets/schema_input.json`](../assets/schema_input.json). See [Usage → Samplesheet input](usage.md#samplesheet-input).
