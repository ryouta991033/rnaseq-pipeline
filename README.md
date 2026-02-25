## RNA-seq Aging Pipeline

A reproducible and modular RNA-seq analysis workflow built with Snakemake.

This pipeline supports:
- Automated SRA download
- Reference genome acquisition (Ensembl)
- Optional read trimming
- Alignment with STAR or HISAT2
- Gene-level quantification
- Differential expression analysis (DESeq2)
- QC aggregation (FastQC + MultiQC)
- Optional AWS S3 result synchronization
Designed for research reproducibility and infrastructure-aware execution.

## Features

- Fully configurable via YAML
- Cloud-ready (AWS S3 integration)
- Modular tool selection (STAR / HISAT2)
- Hybrid input support (SRA + local FASTQ)
- Explicit statistical model control (DESeq2 design & contrast)
- Centralized resource management (thread control)

## Quick Start
1. Clone Repository
```
git clone https://github.com/ryouta991033/rnaseq-pipeline.git
cd rnaseq-aging
```

2. Edit Configuration
Modify:

config/config.yaml
For a full explanation of all parameters, see:

→ docs/configuration.md

3. Run Workflow
```
sudo apt update
↓
sudo apt install -y docker.io
↓
sudo systemctl enable docker 
↓
sudo systemctl start docker 
↓
sudo usermod -aG docker ubuntu 
(If docker is not listed in the output of groups, please run exit once and log in again.)
↓
docker build -t rnaseq-pipeline .
↓
./scripts/setup_env.sh
↓
./scripts/download_fastq.sh
↓
./scripts/generate_samples_tsv.sh
↓
docker run --rm \
-u $(id -u):$(id -g) \
-e HOME=/data \
-v $(pwd):/data \
-w /data \
rnaseq-pipeline \
snakemake --cores 8

```

## Directory Structure

```
rnaseq-pipeline/
│
├── scripts/
│   ├── S3_backup.sh
│   ├── check_s3_result.sh
│   ├── generate_samples_tsv.sh
│   ├── setup_env.sh
│   ├── check_s3_data.sh
│   ├── download_fastq.sh
│   └── run_deseq2.R
│
├── config/
│   └── config.yaml
│
├── data/
│   ├── fastq/
│   └── reference/
│       ├── fasta/
│       └── gtf/
│
├── results/
│   ├── reference/
│   │   ├── STAR_index/
│   │   └── hisat2_index/
│   ├── trim/
│   ├── bam/
│   ├── counts/
│   └── logs/
│
├── workflow/
│   ├── Snakefile
│   └── rules/
│       ├── build_index.smk
│       ├── align.smk
│       ├── trim.smk
│       └── count.smk
│
├── samples.tsv
├── Dockerfile
└── README.md
```

## Workflow Overview
1.SRA download (fasterq-dump)
2.Reference genome download (Ensembl)
3.Optional trimming (Trimmomatic)
4.Alignment (STAR or HISAT2)
5.BAM processing (samtools)
6.Quantification (featureCounts)
7.Differential expression (DESeq2)
8.QC aggregation (MultiQC)
9.Optional S3 synchronization

Example Output

The pipeline produces:
- Alignment BAM files
- Gene count matrix
- DESeq2 results (normalized counts, log2FC, adjusted p-values)
- PCA visualization
- QC reports (HTML)

Example PCA output:
```
results/deseq2/pca_plot.png
```

## Toolchain
Core tools used:
- Snakemake
- STAR / HISAT2
- featureCounts
- DESeq2
- FastQC
- MultiQC
- samtools

## Reference genome source:
- Ensembl release 112

## Reproducibility

- All parameters version-controlled
- Reference versions explicitly defined
- Deterministic DESeq2 contrast configuration
- Centralized thread management

## Optional: AWS Integration
If configured in config.yaml, results can be synchronized to S3.
Ensure AWS CLI is configured before execution.

License
MIT License

Citation
If this pipeline supports your research, please cite appropriately or reference this repository.

## Why This Structure Is Strong

- This README now:
- Is clean and readable
- Does not overwhelm users with config details
- Signals architectural maturity
- Matches structure used in high-quality bioinformatics pipelines
- Looks appropriate for international collaboration
