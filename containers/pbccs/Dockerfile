FROM  mambaorg/micromamba:latest

RUN micromamba create -y -n pbccs-env -c bioconda pbccs && \
    micromamba clean --all --yes

WORKDIR /work

ENTRYPOINT ["micromamba", "run", "-n", "pbccs-env"]