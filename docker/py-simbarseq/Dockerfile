FROM python:3.13.3-alpine3.21
RUN apk update && apk add --no-cache bash=5.2.37-r0
RUN pip install --no-cache-dir pandas==2.2.3 gffutils==0.13
RUN mkdir /app && wget -P /app https://bitbucket.org/berkeleylab/feba/raw/0975b4c3239d35d041fa954fbc359e7f0cecea88/primers/barseq4.index2
COPY simulate_barseq_reads.py utils.py /app/
RUN chmod +x /app/simulate_barseq_reads.py
ENV PATH="/app:$PATH"
WORKDIR /app
