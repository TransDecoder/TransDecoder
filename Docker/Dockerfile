FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

#~~~~~~~~~~~~~~
# Install Preliminaries
#~~~~


RUN apt-get -qq update && apt-get -qq -y install --no-install-recommends \
    automake \
    build-essential \
    bzip2 \
    ca-certificates \
    cpanminus \
    cmake \
    curl \
    default-jre \
    diamond-aligner \
    fort77 \
    ftp \
    g++ \
    gcc \
    gfortran \
    git \
    libblas-dev \
    libbz2-dev \
    libcairo2-dev \
    libcurl4-openssl-dev \
    libdb-dev \
    libghc-zlib-dev \
    libjpeg-dev \
    liblzma-dev \
    libncurses-dev \
    libpcre3-dev \
    libpng-dev \
    libreadline-dev \
    libssl-dev \
    libtbb-dev \
    libx11-dev \
    libxml2-dev \
    libxt-dev \
    libzmq3-dev \
    make \
    ncbi-blast+ \
    nano \
    perl \
    pkg-config \
    python3 \
    python3-dev \
    python3-pip \
    python3-setuptools \
    r-base \
    r-bioc-seqlogo \
    r-cran-caret \
    r-cran-ggplot2 \
    rsync \
    samtools \
    hmmer \
    texlive-latex-base \
    tzdata \
    unzip \
    wget \
    x11-common \
    zlib1g-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*



## Perl stuff
RUN cpanm install URI::Escape

## R packages not available from Ubuntu packages
RUN Rscript -e "install.packages('rpart.plot', repos='https://cloud.r-project.org')"

#~~~~~~~~~~~
# Environmental Variables 
#~~~~~~~~~~~
ENV SRC=/usr/local/src
ENV BIN=/usr/local/bin


# Core bioinformatics tools are installed from Ubuntu packages:
# samtools, ncbi-blast+, diamond-aligner, and hmmer.


# Transdecoder v6.0.0
WORKDIR $SRC
ENV TRANSDECODER_CO=ec18d65529ddacab4cd6c930f8b86406a124b5ea

RUN git clone https://github.com/TransDecoder/TransDecoder.git && \
   cd TransDecoder && \
   git checkout $TRANSDECODER_CO && \
   mv * $BIN
