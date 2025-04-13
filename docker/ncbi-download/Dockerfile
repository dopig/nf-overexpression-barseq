# Stage 1: Builder with Micromamba
FROM mambaorg/micromamba:2.0.8-debian12-slim AS builder
USER root
RUN micromamba create -y --copy -p /env \
    -c conda-forge wget=1.21.4 coreutils=9.5 \
    -c bioconda perl-time-hires=1.9764 entrez-direct=22.4 && \
    micromamba clean --all --yes

# Stage 2: Slim final image
FROM debian:bookworm-20250317-slim AS final
COPY --from=builder /env /env
COPY ncbi_download.sh /app/
RUN apt-get update && apt-get install -y procps=2:4.0.2-3 --no-install-recommends && \
    rm -rf /var/lib/apt/lists/*
RUN chmod +x /app/ncbi_download.sh
ENV PATH="/app:/env/bin:$PATH"
WORKDIR /app
