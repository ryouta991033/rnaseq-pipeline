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
0. As a first step, please attach an IAM role that includes S3AllFreeAccess permissions to the EC2 instance.
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

## Sample Analysis: L4 vs D6 Comparison
Below are representative outputs generated automatically by this RNA-seq pipeline using a test dataset comparing L4 and D6 conditions.

All figures are reproducible and generated directly from raw FASTQ files using the standardized workflow.
PCA Plot

## The PCA plot summarizes global transcriptomic variation across samples.
- Samples cluster according to biological condition (L4 vs D6).
- Replicates group tightly together, indicating good data quality.
- The primary principal component captures condition-specific expression differences.
This confirms that the experimental groups are clearly separated at the transcriptome-wide level.
<img width="1800" height="1500" alt="PCA_plot" src="https://github.com/user-attachments/assets/69692beb-dd4b-4c2c-8a3d-195936f50f40" />

## Heatmap (Top Differentially Expressed Genes)
The heatmap shows the top differentially expressed genes ranked by adjusted p-value.
- Expression values are variance-stabilized (DESeq2 VST).
- Genes are scaled across samples.
- Clear condition-specific expression patterns are observed.
This visualization highlights robust transcriptional differences between L4 and D6.
<img width="1200" height="1600" alt="214f2658-1" src="https://github.com/user-attachments/assets/cd91cdfb-35b3-485c-b7cf-64607d17b584" />

## Volcano Plot
The volcano plot displays:
- log2 Fold Change (x-axis)
- –log10 Adjusted p-value (y-axis)
- Significantly regulated genes are automatically highlighted.
This allows rapid identification of biologically meaningful upregulated and downregulated genes.
<img width="1200" height="1000" alt="34cc059e-1" src="https://github.com/user-attachments/assets/aa0d239f-3884-43eb-aa9c-b2c1f12c501e" />


## Reproducibility
All visualizations are:
- Fully automated
- Docker-compatible
- Generated via Snakemake workflow
- Reproducible from raw FASTQ files
The pipeline is designed to provide publication-ready figures with minimal manual intervention.

This workflow can be customized for external datasets and adapted to specific experimental designs (e.g., multi-group comparisons, time-series, batch correction).
