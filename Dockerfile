FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
       build-essential \
       ca-certificates \
       libboost-all-dev \
       libgmp-dev \
       libmpfi-dev \
       libmpfr-dev \
       python3 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /artifact
COPY . .

RUN sha256sum -c RELEASE_SHA256SUMS \
    && cd proof \
    && sha256sum -c checksums.sha256

CMD ["bash", "scripts/check_release.sh"]
