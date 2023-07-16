FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

#~~~~~~~~~~~~~~
# Install Preliminaries
#~~~~


RUN apt-get -qq update && apt-get -qq -y install \
    automake \
    build-essential \
    bzip2 \
    cmake \
    curl \
    default-jre \
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
    libncurses5-dev \
    libpcre3-dev \
    libpng-dev \
    libreadline-dev \
    libreadline-dev \
    libssl-dev \
    libtbb-dev \
    libx11-dev \
    libxml2-dev \
    libxt-dev \
    libzmq3-dev \
    make \
    nano \
    perl \
    pkg-config \
    python3 \
    python3-dev \
    python3-distutils \
    python3-pip \
    python3-setuptools \
    rsync \
    texlive-latex-base \
    tzdata \
    unzip \
    wget \
    x11-common \
    zlib1g-dev



## Perl stuff
RUN curl -L https://cpanmin.us | perl - App::cpanminus
RUN cpanm install DB_File
RUN cpanm install URI::Escape

#~~~~~~~~~~~
# Environmental Variables 
#~~~~~~~~~~~
ENV SRC /usr/local/src
ENV BIN /usr/local/bin


# samtools
WORKDIR $SRC
ENV SAMTOOLS_VERSION 1.13
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar xvf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && \
    ./configure && make && make install



# blast
WORKDIR $SRC
ENV BLASTPLUS_VERSION 2.12.0
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLASTPLUS_VERSION}/ncbi-blast-${BLASTPLUS_VERSION}+-x64-linux.tar.gz && \
    tar xvf ncbi-blast-${BLASTPLUS_VERSION}+-x64-linux.tar.gz && \
    cp ncbi-blast-${BLASTPLUS_VERSION}+/bin/* $BIN && \
    rm -r ncbi-blast-${BLASTPLUS_VERSION}+



# Hmmer

WORKDIR $SRC
RUN wget http://eddylab.org/software/hmmer/hmmer.tar.gz && \
    tar xvf hmmer.tar.gz && \
    cd hmmer-3.3.2 && \
    ./configure && \
    make && make install





# Trinotate
WORKDIR $SRC
ENV TRANSDECODER_CO=817d772325c8a354aec6b94ac58993b994ac7392

RUN git clone https://github.com/TransDecoder/TransDecoder.git && \
   cd TransDecoder && \
   git checkout $TRANSDECODER_CO && \
   mv * $BIN



