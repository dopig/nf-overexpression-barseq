# Stage 1: Builder with Micromamba
FROM mambaorg/micromamba:2.0.8-debian12-slim AS builder
USER root
COPY environment.yml /opt/env/environment.yml
RUN micromamba create -y --copy -p /env -f /opt/env/environment.yml && \
    micromamba clean --all --yes

# Stage 2: Slim final image
FROM debian:bookworm-20250317-slim AS final
COPY --from=builder /env /env
COPY simulate_barseq_reads.py utils.py /app/
ENV PATH="/opt/bin:/env/bin:$PATH"
WORKDIR /app
