# Start from Rocker R 4.3.3 base image
FROM rocker/r-ver:4.3.3

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        wget \
        curl \
        ca-certificates \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        perl \
        git \
        autoconf \
        automake \
        libtool \
        make \
        gcc \
        g++ \
        bzip2 \
        libncurses5-dev \
        libncursesw5-dev \
        libz-dev \
        && rm -rf /var/lib/apt/lists/*

# Install htslib 1.19 from source
RUN cd /tmp && \
    wget https://github.com/samtools/htslib/releases/download/1.19/htslib-1.19.tar.bz2 && \
    tar -xjf htslib-1.19.tar.bz2 && \
    cd htslib-1.19 && \
    ./configure && \
    make && \
    make install && \
    cd /tmp && rm -rf htslib-1.19*

# Install bedtools from source (latest release)
RUN cd /tmp && \
    git clone https://github.com/arq5x/bedtools2.git && \
    cd bedtools2 && \
    make && \
    cp -r bin/* /usr/local/bin/ && \
    cd /tmp && rm -rf bedtools2

# Set default shell to bash
SHELL ["/bin/bash", "-c"]

# Set working directory
WORKDIR /workspace
