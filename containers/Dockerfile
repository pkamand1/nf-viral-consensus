# Single image containing every tool the pipeline calls.
# Build:  docker build -f containers/Dockerfile -t nf-viral-consensus:0.1.0 .
# Use:    nextflow run main.nf -profile docker --single_container nf-viral-consensus:0.1.0 ...
FROM mambaorg/micromamba:2.0.5

LABEL org.opencontainers.image.title="nf-viral-consensus" \
      org.opencontainers.image.description="Reference-based viral consensus genomes, variant calls, and phylogeny" \
      org.opencontainers.image.authors="Patricia K. Kamanda" \
      org.opencontainers.image.version="0.1.0"

COPY --chown=$MAMBA_USER:$MAMBA_USER envs/environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:$PATH
