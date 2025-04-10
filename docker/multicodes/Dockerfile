# Stage 1: Builder to download files
FROM alpine:3.21.3 AS builder
RUN apk update && apk add --no-cache wget

RUN wget -P /app/src https://bitbucket.org/berkeleylab/feba/raw/0975b4c3239d35d041fa954fbc359e7f0cecea88/bin/MultiCodes.pl
RUN wget -P /app/lib https://bitbucket.org/berkeleylab/feba/raw/0975b4c3239d35d041fa954fbc359e7f0cecea88/lib/FEBA_Utils.pm
RUN wget -P /app/primers https://bitbucket.org/berkeleylab/feba/raw/0975b4c3239d35d041fa954fbc359e7f0cecea88/primers/barseq4.index2

# Stage 2: Final image with necessary runtime dependencies
FROM alpine:3.21.3
RUN apk update && apk add --no-cache bash=5.2.37-r0 perl=5.40.1-r0

# Copy the downloaded files from the builder stage
COPY --from=builder /app /app
RUN chmod +x /app/src/MultiCodes.pl