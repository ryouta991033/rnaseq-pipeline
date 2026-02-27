# =====================================================
# Base image
# =====================================================
FROM ubuntu:22.04

LABEL maintainer="User"
LABEL description="RNA-seq pipeline container (Snakemake + STAR + HISAT2 + DESeq2)"
LABEL version="1.0"

ENV DEBIAN_FRONTEND=noninteractive
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

# =====================================================
# System packages
# =====================================================
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    wget \
    git \
    unzip \
    zip \
    gzip \
    pigz \
    parallel \
    bedtools \
    htop \
    vim \
    less \
    openjdk-17-jdk \
    python3 \
    python3-pip \
    bzip2 \
    ca-certificates \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# =====================================================
# Miniconda
# =====================================================
ENV CONDA_DIR=/opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH

RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p $CONDA_DIR && \
    rm /tmp/miniconda.sh

# Conda TOS
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Strict channel priority（重要）
RUN conda config --set channel_priority strict

# =====================================================
# Core bioinformatics tools (conda統一)
# =====================================================
RUN conda install -y -c conda-forge -c bioconda \
    snakemake \
    star \
    hisat2 \
    samtools \
    subread \
    stringtie \
    multiqc \
    yq \
    parallel \
    pigz \
    r-base \
    bioconductor-deseq2 \
    bioconductor-tximport \
    r-optparse \
    r-tidyverse \
    r-pheatmap \
    r-ggplot2 \
    r-ggrepel \
    bioconductor-rtracklayer \
 && conda clean -afy

# =====================================================
# Python libraries
# =====================================================
RUN pip install --no-cache-dir numpy pandas

# =====================================================
# AWS CLI v2
# =====================================================
RUN curl -s https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip -o /tmp/awscliv2.zip && \
    unzip -q /tmp/awscliv2.zip -d /tmp && \
    /tmp/aws/install && \
    rm -rf /tmp/aws /tmp/awscliv2.zip

# =====================================================
# Trimmomatic
# =====================================================
WORKDIR /opt
RUN wget -q http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip \
 && unzip Trimmomatic-0.39.zip \
 && ln -s /opt/Trimmomatic-0.39/trimmomatic-0.39.jar /usr/local/bin/trimmomatic.jar \
 && rm Trimmomatic-0.39.zip

RUN echo '#!/bin/bash\nexec java -jar /usr/local/bin/trimmomatic.jar "$@"' \
 > /usr/local/bin/trimmomatic \
 && chmod +x /usr/local/bin/trimmomatic

# =====================================================
# ncbi-sra-tools
# =====================================================
RUN wget -q https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.3/sratoolkit.2.11.3-ubuntu64.tar.gz \
 && tar -xzf sratoolkit.2.11.3-ubuntu64.tar.gz -C /opt \
 && rm sratoolkit.2.11.3-ubuntu64.tar.gz

ENV PATH="/opt/sratoolkit.2.11.3-ubuntu64/bin:${PATH}"


# =====================================================
# Default workdir
# =====================================================
WORKDIR /data

CMD ["bash"]

