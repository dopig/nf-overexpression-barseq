# Stage 1: Builder with Micromamba
FROM mambaorg/micromamba:2.0.8-debian12-slim AS builder
USER root

RUN micromamba create -y --copy -p /env \
    -c conda-forge python=3.13.3 pandas==2.2.3 matplotlib==3.10.1 \
    -c bioconda minimap2=2.29 vsearch=2.30.0 \
    && micromamba clean --all --yes

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget=1.21.3-1+deb12u1 \
    git=1:2.39.5-0+deb12u2 \
    && rm -rf /var/lib/apt/lists/*

# Clone Boba-seq repository and fix issues with current pandas, current matplotlib, and to avoid docstring syntax errors
RUN git clone --depth 1 https://github.com/OGalOz/Boba-seq.git /app && \
    sed -i \
        -e 's/pos_bc = pos_bc.append(df_bc_keep\[\["query", "qhi_x", "qlo_y"\]\])/pos_bc = pd.concat(\[pos_bc, df_bc_keep\[\["query", "qhi_x", "qlo_y"\]\]\])/g' \
        -e 's/pos_ins = pos_ins.append(df_ins_keep\[\["query", "qhi_x", "qlo_y"\]\])/pos_ins = pd.concat(\[pos_ins, df_ins_keep\[\["query", "qhi_x", "qlo_y"\]\]\])/g' \
        /app/src/step2.py && \
    sed -i 's/"seaborn-whitegrid"/"seaborn-v0_8-whitegrid"/g' /app/src/step6.py && \
    sed -i -e '1 s/"""/r"""/' -e '130 s/"""/r"""/' -e '309 s/"""/r"""/' /app/src/step1.py && \
    sed -i -e '1 s/"""/r"""/' -e '108 s/"""/r"""/' /app/src/step3.py

# Install usearch (adjust version as needed)
RUN wget https://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz && \
    gunzip usearch11.0.667_i86linux32.gz && \
    chmod +x usearch11.0.667_i86linux32 && \
    mv usearch11.0.667_i86linux32 /app/usearch


# Stage 2: Slim final image
FROM debian:bookworm-20250317-slim AS final

RUN apt-get update && apt-get install -y --no-install-recommends \
    procps=2:4.0.2-3 && \
    rm -rf /var/lib/apt/lists/*

COPY --from=builder /env /env
COPY --from=builder /app /app
ENV PATH="/app:/env/bin:$PATH"

COPY ./adapt_json.py /opt/bin/adapt_json.py

RUN mkdir -p /data/fastq

WORKDIR /work
