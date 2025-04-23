# Stage 1: Builder with Micromamba
FROM debian:bookworm-20250317-slim AS builder
USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget=1.21.3-1+deb12u1 ca-certificates=20230311 && \
    rm -rf /var/lib/apt/lists/*
RUN wget -P /app/src https://bitbucket.org/berkeleylab/feba/raw/0975b4c3239d35d041fa954fbc359e7f0cecea88/bin/MultiCodes.pl
RUN wget -P /app/lib https://bitbucket.org/berkeleylab/feba/raw/0975b4c3239d35d041fa954fbc359e7f0cecea88/lib/FEBA_Utils.pm
RUN wget -P /app/primers https://bitbucket.org/berkeleylab/feba/raw/0975b4c3239d35d041fa954fbc359e7f0cecea88/primers/barseq4.index2


# Stage 2: Slim final image
FROM debian:bookworm-20250317-slim AS final
RUN apt-get update && apt-get install -y --no-install-recommends procps=2:4.0.2-3 perl=5.36.0-7+deb12u2 \
    && rm -rf /var/lib/apt/lists/*
COPY --from=builder /app /app
RUN chmod +x /app/src/MultiCodes.pl
