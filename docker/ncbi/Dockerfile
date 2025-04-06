FROM  mambaorg/micromamba:latest

RUN micromamba create -y -n ncbi-env -c conda-forge wget coreutils -c bioconda perl-time-hires entrez-direct && \
    micromamba clean --all --yes

WORKDIR /work

ENTRYPOINT ["micromamba", "run", "-n", "ncbi-env"]