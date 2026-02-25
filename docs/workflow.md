Workflow Overview

This document describes the logical structure and execution flow of the RNA-seq pipeline.

The workflow is implemented using Snakemake and is fully configuration-driven.
All behavior is controlled via config/config.yaml.

Pipeline Execution Graph
The workflow follows a linear but modular structure:
```
SRA Download / User FASTQ
        ↓
Reference Genome Acquisition
        ↓
(Optional) Read Trimming
        ↓
Alignment (STAR or HISAT2)
        ↓
BAM Processing (samtools)
        ↓
Quantification (featureCounts)
        ↓
Differential Expression (DESeq2)
        ↓
(Optional) AWS S3 Sync
```
Each stage is conditionally executed based on configuration parameters.

Step-by-Step Description

1. Data Acquisition
Public Data
- FASTQ files are downloaded from SRA using fasterq-dump.
Local Data
- If userdata is defined, local FASTQ files are used instead.
This design enables hybrid workflows combining public and private datasets.

2. Reference Genome Setup

If not already present, the pipeline:
- Downloads reference FASTA and GTF from Ensembl
- Builds index for selected aligner (STAR or HISAT2)
Reference versioning is explicitly controlled via configuration.

3. Optional Read Trimming

If enabled (trim.enable: true):
- Adapter trimming and quality filtering are performed using Trimmomatic.
- Output FASTQ replaces raw reads for downstream alignment.
If disabled, raw FASTQ files are used directly.

4. Alignment
The mapper is selected via:
```
mapper: star  # or hisat2
```
STAR
- Genome index built automatically if absent
- Produces coordinate-sorted BAM

HISAT2
- Index built automatically
- BAM sorted using samtools
The abstraction layer ensures mapper switching without altering workflow code.

5. BAM Processing

- Sorting (if required)
- Index generation
- Storage under results/bam/
This stage prepares files for quantification.

6. Quantification

Gene-level counts are generated using featureCounts.
Configurable parameters:
- Strandedness
- Paired-end mode
- Additional featureCounts options
Output:
```
results/counts/gene_counts.tsv
```

7. Differential Expression Analysis

DE analysis is performed using DESeq2.
Configuration-driven:
```
deseq2:
  design: "~ condition"
  contrast:
    factor: "condition"
    numerator: "L4"
    denominator: "D6"
```

Outputs:
- Normalized counts
- Log2 fold changes
- Adjusted p-values
- PCA plot
- MA plot
Stored in:
```
results/deseq2/
```

8. Quality Control

If enabled:
- FastQC runs per sample
- MultiQC aggregates reports
Outputs stored in:
```
results/qc/
results/multiqc/
```

9. Optional Cloud Synchronization

If AWS parameters are configured:
- Results are synchronized to S3
- Controlled by terminal rule
This allows cloud-backed reproducibility and sharing.

Modular Design Principles

The workflow is structured around:
- Parameter-driven execution
- Conditional rule activation
- Tool abstraction layer
- Explicit input/output declaration
- Deterministic statistical modeling

No modification to the Snakefile is required when:
- Adding samples
- Changing aligners
- Updating reference genome
- Modifying statistical design
- Scaling computational resources

Reproducibility Strategy

Reproducibility is ensured by:
- Version-controlled configuration
- Explicit reference genome URLs
- Centralized thread management
- Deterministic DESeq2 contrast specification
- Structured output directories

Extensibility

The pipeline can be extended by:
- Adding new aligners
- Integrating transcript-level quantification
- Adding alternative DE methods
- Implementing cluster profiles (SLURM / SGE)
- Expanding cloud deployment strategies
The current architecture supports modular expansion.

Execution Example
```
snakemake --cores 8
```
With Conda:
```
snakemake --use-conda --cores 8
```
Inside Docker:
```
docker compose up
```

Summary

This pipeline is designed not only as an RNA-seq workflow, but as a configurable, reproducible analytical framework suitable for:
Academic research
- Collaborative projects
- Cloud-enabled environments
- Infrastructure-aware bioinformatics workflows
