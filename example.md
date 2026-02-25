Reproducible RNA-seq Analysis Pipeline

Example analysis: Aging transition (L4 vs Day 6) in Caenorhabditis elegans

This repository demonstrates a fully reproducible bulk RNA-seq workflow using Snakemake and Docker.
The pipeline performs quality control, alignment, quantification, differential expression analysis (DESeq2), and visualization.

🔬 Example Dataset

Organism: Caenorhabditis elegans
Comparison: L4 vs Day 6 (early aging transition)
Replicates: 2 vs 2
Data source: NCBI SRA

📊 Key Results
Principal Component Analysis (PCA)

Clear separation between L4 and Day 6 samples
Replicates cluster together
Indicates global transcriptomic remodeling during early aging
<img width="1800" height="1500" alt="DESeq2_PCA" src="https://github.com/user-attachments/assets/7c807d5f-127f-4cf0-b9a0-a20aa4cca92e" />
