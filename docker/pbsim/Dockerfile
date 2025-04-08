# Stage 1: Builder with Micromamba
FROM mambaorg/micromamba:2.0.8-debian12-slim AS builder
USER root
RUN micromamba create -y --copy -p /env \
    -c conda-forge -c bioconda pbsim3=3.0.5 samtools=1.21 && \
    micromamba clean --all --yes

# Stage 2: Slim final image
FROM debian:bookworm-20250317-slim AS final
COPY --from=builder /env /env
ENV PATH="/env/bin:$PATH"
WORKDIR /work
