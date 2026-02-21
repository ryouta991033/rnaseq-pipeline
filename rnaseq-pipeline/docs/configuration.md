## Configuration Guide

All runtime behavior of the RNA-seq pipeline is controlled via:
```
config/config.yaml
```
The workflow logic is intentionally separated from configuration to ensure:
- Reproducibility
- Transparency
- Infrastructure portability
- Experimental flexibility without modifying workflow code

Configuration Structure Overview
```
project
aws
download
samples
userdata (optional)
mapper
qc
workdir / outdir
threads
trim
raw
reference
featurecounts
deseq2
results
logs
```
1. Project Metadata
project.name

Project identifier used for labeling outputs and cloud prefixes.
```
project:
  name: rnaseq-aging
```
2. AWS Configuration (Optional)

Used only if S3 synchronization is enabled.
```
aws:
  bucket: bioinfo-project
  region: ap-northeast-1
  results_prefix: rnaseq-aging
```

| Parameter        | Description                            |
| ---------------- | -------------------------------------- |
| `bucket`         | S3 bucket name                         |
| `region`         | AWS region                             |
| `results_prefix` | Project subdirectory inside the bucket |

3. Data Download
download.raw_fastq_dir
Directory where downloaded FASTQ files are stored.

Reference Genome Download
```
download:
  reference:
    base_dir:
    fasta_dir:
    gtf_dir:
```
FASTA
```
fasta:
  url:
  name:
```
GTF
```
gtf:
  url:
  name:
```
This design ensures:
- Versioned reference tracking
- Automatic reproducibility
- No hidden manual dependency

4. Sample Definition

Each sample is defined as:
```
samples:
  sample1:
    sra: SRRxxxxxxx
    condition: L4
    layout: SE
```
| Field       | Description            |
| ----------- | ---------------------- |
| `sra`       | SRA accession ID       |
| `condition` | Experimental condition |
| `layout`    | `SE` or `PE`           |

The condition field is used for DESeq2 modeling.

5. User-Provided FASTQ (Optional)

To include local FASTQ files instead of SRA:
```
userdata:
  user_sample1:
    layout: SE
    fq1: data/fastq/sample.fastq.gz
    condition: D6
```
Supports hybrid workflows (public + private datasets).

6. Mapper Selection
```
mapper: star   # star | hisat2
```

Supported options:
- star
- hisat2
Index building is handled automatically.
7. Quality Control
```
qc:
  fastqc: true
  multiqc: true
```
Boolean flags to enable/disable QC steps.

8. Working Directories
```
workdir: data
outdir: results
```
Allows flexible data/result relocation.

9. Thread Management

Centralized resource allocation:
```
threads:
  fasterq: 4
  pigz: 8
  trimmomatic: 4
  star: 8
  hisat2: 8
  samtools: 4
  featurecounts: 4
```
Improves scalability and cluster compatibility.

10. Trimming
```
trim:
  enable: true
  trimmomatic_opts: "SLIDINGWINDOW:4:20 MINLEN:36"
  adapters: "adapters/TruSeq3-PE.fa"
```
| Parameter          | Description              |
| ------------------ | ------------------------ |
| `enable`           | Toggle trimming          |
| `trimmomatic_opts` | Quality trimming options |
| `adapters`         | Adapter file path        |

11. Raw & Reference Paths

These fields explicitly declare runtime file locations.

Raw
```
raw:
  fastq:
  metadata:
```
Reference
```
reference:
  fasta:
  gtf:
  star_index:
  hisat2_index:
  read_length:
```
Explicit declaration prevents hidden path assumptions.

12. Quantification (featureCounts)
```
featurecounts:
  stranded: 0
  paired: false
  extra_opts: ""
```
| Value | Meaning          |
| ----- | ---------------- |
| 0     | Unstranded       |
| 1     | Stranded         |
| 2     | Reverse-stranded |

13. Differential Expression (DESeq2)
```
deseq2:
  design: "~ condition"
  contrast:
    factor: "condition"
    numerator: "L4"
    denominator: "D6"
  reference:
    condition: "L4"
```
design

Statistical formula passed to DESeq2.

contrast

Defines comparison:
- numerator = treatment
- denominator = control

reference
Factor releveling for deterministic modeling.

14. Result Directories
```
results:
  base:
  bam:
  counts:
  deseq2:
  qc:
  multiqc:
  logs:
```
All outputs are centralized and configurable.

15. Logging
```
logs:
  snakemake:
```
Stores execution logs for reproducibility and debugging.

Design Philosophy

This configuration system emphasizes:
- Separation of workflow and parameters
- Deterministic statistical modeling
- Infrastructure awareness
- Explicit resource control
- Version-traceable reference data

It is structured for:
- Academic reproducibility
- Collaborative projects
- Cloud execution
- Production-like deployment scenarios

Best Practices

- Keep config.yaml version-controlled
- Use config.example.yaml for template sharing
- Avoid modifying the Snakefile directly
- Track reference genome versions explicitly

Summary

The pipeline is fully driven by configuration.
No workflow modification is required when:
- Changing samples
- Switching aligners
- Updating reference genome
- Modifying statistical contrast
- Scaling compute resources


