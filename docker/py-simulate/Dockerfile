# Stage 1: Builder with Micromamba
FROM mambaorg/micromamba:2.0.8-debian12-slim AS builder
USER root
RUN micromamba create -y --copy -p /env \
    -c conda-forge python=3.13.3 numpy==2.2.3 scipy==1.15.2 pandas==2.2.3 \
    -c bioconda biopython==1.85 gffutils==0.13 \
    && micromamba clean --all --yes
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget=1.21.3-1+deb12u1 \
    && rm -rf /var/lib/apt/lists/*
RUN wget -P /app https://bitbucket.org/berkeleylab/feba/raw/0975b4c3239d35d041fa954fbc359e7f0cecea88/primers/barseq4.index2

# Stage 2: Slim final image
FROM debian:bookworm-20250317-slim AS final
RUN apt-get update && apt-get install -y --no-install-recommends \
    procps=2:4.0.2-3 && \
    rm -rf /var/lib/apt/lists/*

COPY simulate_library.py simulate_barseq_reads.py utils.py /app/
RUN chmod +x /app/simulate_library.py /app/simulate_barseq_reads.py
COPY --from=builder /app/barseq4.index2 /app/
COPY --from=builder /env /env

ENV PATH="/app:/env/bin:$PATH"
