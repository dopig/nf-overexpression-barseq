FROM python:3.13.3-slim-bookworm
RUN pip install --no-cache-dir numpy==2.2.3 scipy==1.15.2 biopython==1.85
COPY simulate_library.py utils.py /app/
RUN chmod +x /app/simulate_library.py
RUN apt-get update && apt-get install -y procps=2:4.0.2-3 --no-install-recommends && \
    rm -rf /var/lib/apt/lists/*
ENV PATH="/app:$PATH"
WORKDIR /app