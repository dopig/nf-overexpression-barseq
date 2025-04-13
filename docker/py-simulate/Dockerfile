FROM python:3.13.3-slim-bookworm

# System utilities (combine needs of both)
RUN apt-get update && apt-get install -y --no-install-recommends \
    procps=2:4.0.2-3 \
    wget=1.21.3-1+deb12u1 \
    && rm -rf /var/lib/apt/lists/*

# Python packages
RUN pip install --no-cache-dir \
    numpy==2.2.3 \
    scipy==1.15.2 \
    biopython==1.85 \
    pandas==2.2.3 \
    gffutils==0.13

# App directory + shared utils
COPY simulate_library.py simulate_barseq_reads.py utils.py /app/
RUN chmod +x /app/simulate_library.py /app/simulate_barseq_reads.py

# Extra file needed for barseq
RUN wget -P /app https://bitbucket.org/berkeleylab/feba/raw/0975b4c3239d35d041fa954fbc359e7f0cecea88/primers/barseq4.index2

ENV PATH="/app:$PATH"
WORKDIR /app
