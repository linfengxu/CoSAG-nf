FROM mambaorg/micromamba:latest

LABEL maintainer="Linfeng Xu"
LABEL description="CoSAG-nf: CheckM2 1.0.2 + SPAdes 4.2.0"
LABEL version="1.0.0"

USER root

RUN apt-get update && \
    apt-get install -y --no-install-recommends procps && \
    rm -rf /var/lib/apt/lists/*

RUN micromamba create -y -n spades \
    -c bioconda -c conda-forge \
    spades=3.15.5 && \
    micromamba clean --all --yes

RUN micromamba create -y -n checkm2 \
    -c bioconda -c conda-forge \
    checkm2=1.0.2 numpy scipy biopython && \
    micromamba clean --all --yes

ENV PATH="/opt/conda/envs/checkm2/bin:/opt/conda/envs/spades/bin:/usr/local/bin:$PATH"
ENV CONDA_DEFAULT_ENV="checkm2"

RUN checkm2 --version && \
    spades.py --version && \
    python -c "import numpy, scipy, Bio; print('All deps OK')"
