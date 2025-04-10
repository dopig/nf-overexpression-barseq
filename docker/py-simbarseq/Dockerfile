FROM python:3.13.3-alpine3.21
RUN apk update && apk add --no-cache bash=5.2.37-r0
RUN pip install --no-cache-dir pandas==2.2.3 gffutils==0.13
COPY simulate_barseq_reads.py utils.py /app/
ENV PATH="/app:$PATH"
WORKDIR /app