import os

REF_DIR = "results/reference"
THREADS = config.get("threads", {})
FASTA = config["reference"]["fasta"]

################################
# STAR index
################################

STAR_INDEX_DIR = os.path.join(REF_DIR, "STAR_index")

rule build_star_index:
    input:
        fasta=FASTA
    output:
        os.path.join(STAR_INDEX_DIR, "Genome")
    threads: THREADS.get("star", 8)
    shell:
        r"""
        mkdir -p {STAR_INDEX_DIR}
        gunzip -c {input.fasta} > {STAR_INDEX_DIR}/genome.fa
        STAR --runMode genomeGenerate \
             --genomeDir {STAR_INDEX_DIR} \
             --genomeFastaFiles {STAR_INDEX_DIR}/genome.fa \
             --runThreadN {threads}
        """

################################
# HISAT2 index
################################

HISAT2_INDEX_DIR = os.path.join(REF_DIR, "hisat2_index")
HISAT2_INDEX_PREFIX = os.path.join(HISAT2_INDEX_DIR, "genome")

rule build_hisat2_index:
    input:
        fasta=FASTA
    output:
        HISAT2_INDEX_PREFIX + ".1.ht2"
    threads: THREADS.get("hisat2", 8)
    shell:
        r"""
        mkdir -p {HISAT2_INDEX_DIR}
        gunzip -c {input.fasta} > {HISAT2_INDEX_DIR}/genome.fa
        hisat2-build -p {threads} \
            {HISAT2_INDEX_DIR}/genome.fa \
            {HISAT2_INDEX_PREFIX}
        rm {HISAT2_INDEX_DIR}/genome.fa
        """

