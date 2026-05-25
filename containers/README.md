# CoSAG-nf container images

## Published images (Quay.io)

| Image            | Pull                                             | Used for                                            |
| ---------------- | ------------------------------------------------ | --------------------------------------------------- |
| Python bio       | `docker pull quay.io/xulf2022/python3.8_bio:v1`  | Clustering, JSON merge, HTML report (`bin/*.py`)    |
| SPAdes + CheckM2 | `docker pull quay.io/xulf2022/spades_checkm2:v1` | TNF co-assembly optimization (`COSAG_OPTIMIZATION`) |

Image URIs are set in each `modules/local/*/main.nf` `container` line.

## Biocontainers (assembly & databases)

These processes use [Biocontainers](https://biocontainers.pro/) on Quay (pulled automatically with `-profile docker` or `singularity`):

| Tool     | Typical image                                       |
| -------- | --------------------------------------------------- |
| SPAdes   | `quay.io/biocontainers/spades:3.15.5--h95f258a_1`   |
| CheckM2  | `quay.io/biocontainers/checkm2:1.0.2--pyh7cba7a3_0` |
| Sourmash | `quay.io/biocontainers/sourmash:4.8.4--hdfd78af_0`  |
| GTDB-Tk  | `quay.io/biocontainers/gtdbtk:2.4.0--pyhdfd78af_1`  |
| barrnap  | `quay.io/biocontainers/barrnap:0.9--0`              |

Override any image in [`conf/modules.config`](../conf/modules.config).

## Rebuild custom images (optional)

```bash
# From repository root
docker build -t quay.io/xulf2022/python3.8_bio:v1 -f containers/python-bio/Dockerfile containers/python-bio
docker build -t quay.io/xulf2022/spades_checkm2:v1 -f containers/spades-checkm2/Dockerfile containers/spades-checkm2

docker push quay.io/xulf2022/python3.8_bio:v1
docker push quay.io/xulf2022/spades_checkm2:v1
```

## Singularity / Apptainer (HPC)

```bash
singularity pull docker://quay.io/xulf2022/python3.8_bio:v1
singularity pull docker://quay.io/xulf2022/spades_checkm2:v1
```

Then run with `-profile singularity`.

See the main [README](../README.md) for full usage.
