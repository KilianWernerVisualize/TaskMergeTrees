#!/usr/bin/env bash
#Assume Ubuntu Linux 18.04.3
apt-get update \
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
        libvtk7-dev=7.1.1+dfsg1-2 \
        ninja-build=1.8.2-1 \
        python=2.7.15~rc1-1 \
        unzip=6.0-21ubuntu1 \
        zlib1g-dev=1:1.2.11.dfsg-0ubuntu2 \
    && rm -rf /var/lib/apt/lists/*

set -e # abort on error
set -x # verbose output

# -- Teem

export TEEM_VERSION="1.11.0"

mkdir teem-source

curl -L https://netcologne.dl.sourceforge.net/project/teem/teem/1.11.0/teem-1.11.0-src.tar.gz | tar zx -C teem-source --strip-components 1

mkdir teem-build
pushd teem-build

cmake -G Ninja \
      -D CMAKE_INSTALL_PREFIX=/usr \
      ../teem-source

cmake --build . --target install

popd

# rm -rf teem-source teem-build

# --- HPX

export HPX_VERSION="1.4.1"

mkdir hpx-source

curl -L https://github.com/STEllAR-GROUP/hpx/archive/${HPX_VERSION}.tar.gz | tar zx -C hpx-source --strip-components 1

mkdir hpx-build
pushd hpx-build

cmake -G Ninja \
      -D CMAKE_INSTALL_PREFIX=/usr \
      -D CMAKE_BUILD_TYPE=Release \
      -D HPX_WITH_MALLOC=tcmalloc \
      -D HPX_WITH_EXAMPLES=OFF \
      -D HPX_WITH_TESTS=OFF \
      -D HPX_WITH_PAPI=ON \
      ../hpx-source

cmake --build . --target install

popd

rm -rf hpx-source hpx-build

set -ex

mkdir build
pushd build

cmake -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      code/src

cmake --build .

./hpxct \
    /data/ctBones.vti \
   	| tee /results/output

ls
cp output_0.vti /results/
cp output.vtp /results/

popd
