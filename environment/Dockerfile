# hash:sha256:f555575891b49e3d8f20397cf3423c5ed4b1fb702d6f8bbf808afcce8c022d5f
FROM registry.codeocean.com/codeocean/ubuntu:18.04.3

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential=12.4ubuntu1 \
        cmake=3.10.2-1ubuntu2.18.04.1 \
        curl=7.58.0-2ubuntu3.9 \
        git=1:2.17.1-1ubuntu0.7 \
        libblosc-dev=1.14.2+ds1-1 \
        libboost-atomic-dev=1.65.1.0ubuntu1 \
        libboost-chrono-dev=1.65.1.0ubuntu1 \
        libboost-date-time-dev=1.65.1.0ubuntu1 \
        libboost-filesystem-dev=1.65.1.0ubuntu1 \
        libboost-program-options-dev=1.65.1.0ubuntu1 \
        libboost-regex-dev=1.65.1.0ubuntu1 \
        libboost-system-dev=1.65.1.0ubuntu1 \
        libboost-thread-dev=1.65.1.0ubuntu1 \
        libbz2-dev=1.0.6-8.1ubuntu0.2 \
        libglm-dev=0.9.9~a2-2 \
        libgoogle-perftools-dev=2.5-2.2ubuntu3 \
        libhwloc-dev=1.11.9-1 \
        libopenmpi-dev=2.1.1-8 \
        libpapi-dev=5.6.0-1 \
        libpng-dev=1.6.34-1ubuntu0.18.04.2 \
        libsnappy-dev=1.1.7-1 \
        libtbb-dev=2017~U7-8 \
        libvtk7-dev \
        ninja-build=1.8.2-1 \
        python=2.7.15~rc1-1 \
        unzip=6.0-21ubuntu1 \
        zlib1g-dev=1:1.2.11.dfsg-0ubuntu2 \
    && rm -rf /var/lib/apt/lists/*

COPY postInstall /
RUN /postInstall
