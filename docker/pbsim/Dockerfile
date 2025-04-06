FROM  mambaorg/micromamba:latest

RUN micromamba create -y -n pbsim-env -c conda-forge -c bioconda pbsim3 samtools && \
    micromamba clean --all --yes

WORKDIR /work

ENTRYPOINT ["micromamba", "run", "-n", "pbsim-env"]