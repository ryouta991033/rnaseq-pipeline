################################
# align.smk（index依存対応版）
################################

import os

RESULTS_DIR = "results"
MAP_DIR = os.path.join(RESULTS_DIR, "bam")
TRIM_DIR = os.path.join(RESULTS_DIR, "trim")
REF_DIR = os.path.join(RESULTS_DIR, "reference")

MAPPER = config.get("mapper", "star")

THREADS = {
    "star": config.get("threads", {}).get("star", 8),
    "hisat2": config.get("threads", {}).get("hisat2", 8),
}

################################
# STAR
################################

if MAPPER == "star":

    STAR_INDEX_DIR = os.path.join(REF_DIR, "STAR_index")

    rule align:
        input:
            index=os.path.join(STAR_INDEX_DIR, "Genome"),
            r1=os.path.join(TRIM_DIR, "{sample}.trim.fastq.gz")
        output:
            bam=os.path.join(MAP_DIR, "{sample}.bam"),
            bai=os.path.join(MAP_DIR, "{sample}.bam.bai")
        threads: THREADS["star"]
        shell:
            r"""
            mkdir -p {MAP_DIR}

            STAR \
                --genomeDir {STAR_INDEX_DIR} \
                --readFilesIn {input.r1} \
                --readFilesCommand zcat \
                --runThreadN {threads} \
                --outSAMtype BAM SortedByCoordinate \
                --outFileNamePrefix {MAP_DIR}/{wildcards.sample}.

            mv {MAP_DIR}/{wildcards.sample}.Aligned.sortedByCoord.out.bam {output.bam}

            samtools index {output.bam}
            """

################################
# HISAT2
################################

elif MAPPER == "hisat2":

    HISAT2_INDEX_PREFIX = os.path.join(REF_DIR, "hisat2_index", "genome")

    rule align:
        input:
            index=HISAT2_INDEX_PREFIX + ".1.ht2",
            r1=os.path.join(TRIM_DIR, "{sample}.trim.fastq.gz")
        output:
            bam=os.path.join(MAP_DIR, "{sample}.bam"),
            bai=os.path.join(MAP_DIR, "{sample}.bam.bai")
        threads: THREADS["hisat2"]
        shell:
            r"""
            mkdir -p {MAP_DIR}

            hisat2 \
                -x {HISAT2_INDEX_PREFIX} \
                -U {input.r1} \
                -p {threads} \
            | samtools sort -@ {threads} -o {output.bam}

            samtools index {output.bam}
            """

