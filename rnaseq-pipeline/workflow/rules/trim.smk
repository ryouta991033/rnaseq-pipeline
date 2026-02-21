################################
# trim.smk（最終安定・統一版）
################################

import os
import pandas as pd

TRIM_DIR = "results/trim"

samples_df = pd.read_csv("samples.tsv", sep="\t").set_index("sample")

# -----------------------------
# layout取得
# -----------------------------
def get_layout(wc):
    return samples_df.loc[wc.sample, "layout"]

def get_fq1(wc):
    return str(samples_df.loc[wc.sample, "fq1"])

def get_fq2(wc):
    fq2 = samples_df.loc[wc.sample, "fq2"]
    return str(fq2) if fq2 != "-" else None


################################
# PE
################################
rule trimmomatic_pe:
    input:
        r1=get_fq1,
        r2=get_fq2
    output:
        r1=os.path.join(TRIM_DIR, "{sample}.R1.trim.fastq.gz"),
        r2=os.path.join(TRIM_DIR, "{sample}.R2.trim.fastq.gz")
    threads: config.get("threads", {}).get("trimmomatic", 4)
    run:
        if get_layout(wildcards) != "PE":
            return

        shell("mkdir -p {TRIM_DIR}")

        shell(r"""
            trimmomatic PE \
                -threads {threads} \
                {input.r1} {input.r2} \
                {output.r1} {TRIM_DIR}/{wildcards.sample}.R1.unpaired.fastq.gz \
                {output.r2} {TRIM_DIR}/{wildcards.sample}.R2.unpaired.fastq.gz \
                ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10 \
                SLIDINGWINDOW:4:20 \
                MINLEN:36
        """)


################################
# SE
################################
rule trimmomatic_se:
    input:
        r1=get_fq1
    output:
        os.path.join(TRIM_DIR, "{sample}.trim.fastq.gz")
    threads: config.get("threads", {}).get("trimmomatic", 4)
    run:
        if get_layout(wildcards) != "SE":
            return

        shell("mkdir -p {TRIM_DIR}")

        shell(r"""
            trimmomatic SE \
                -threads {threads} \
                {input.r1} \
                {output} \
                ILLUMINACLIP:adapters/TruSeq3-SE.fa:2:30:10 \
                SLIDINGWINDOW:4:20 \
                MINLEN:36
        """)

