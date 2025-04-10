FROM python:3.13.3-slim-bookworm
RUN pip install --no-cache-dir numpy==2.2.3 scipy==1.15.2 biopython==1.85
COPY simulate_library.py utils.py /app/
RUN chmod +x /app/simulate_library.py
ENV PATH="/app:$PATH"
WORKDIR /app