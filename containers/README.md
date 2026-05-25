# CoSAG-nf container images

Published images on [Quay.io](https://quay.io/repository/xulf2022):

| Image | Pull | Used for |
|-------|------|----------|
| Python bio | `docker pull quay.io/xulf2022/python3.8_bio:v1` | Clustering, JSON merge, HTML report (`bin/*.py`) |
| SPAdes + CheckM2 | `docker pull quay.io/xulf2022/spades_checkm2:v1` | TNF co-assembly optimization (`COSAG_OPTIMIZATION`) |

Image URIs are set directly in each `modules/local/*/main.nf` process `container` line.

## Build from Dockerfiles (optional)

Dockerfiles live in this directory if you need to rebuild or retag:

```
containers/
├── python-bio/Dockerfile
├── spades-checkm2/Dockerfile
└── README.md
```

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

Then run with `-profile singularity`; Nextflow will pull or use cached images from the registry.

Other tools (standalone SPAdes, CheckM2 batch, Sourmash, GTDB-Tk, barrnap) still use local `.sif` paths in their modules until additional images are published.

See the main [README](../README.md) for full usage.
