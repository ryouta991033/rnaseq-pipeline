################################
# featureCounts + DESeq2
################################

import os
import pandas as pd

samples = pd.read_csv("samples.tsv", sep="\t")
SAMPLES = samples["sample"].tolist()

# =============================
# featureCounts
# =============================
rule featurecounts:
    input:
        bams=lambda wildcards: [f"results/bam/{s}.bam" for s in SAMPLES],
        gtf=config["reference"]["gtf"]
    output:
        counts="results/counts/gene_counts.txt"
    threads: config["threads"].get("featurecounts", 4)
    shell:
        """
        mkdir -p results/counts
        featureCounts \
            -T {threads} \
            -a {input.gtf} \
            -o {output.counts} \
            -t exon \
            -g gene_id \
            {input.bams}
        """

# =============================
# DESeq2 + gene_name付加
# =============================
rule deseq2:
    input:
        counts="results/counts/gene_counts.txt"
    output:
        all_results="results/deseq2/DESeq2_results_all.csv",
        sig_results="results/deseq2/DESeq2_results_significant.csv",
        rld_object="results/deseq2/rld.rds"
    params:
        conditions=lambda wildcards: samples.set_index("sample") \
                                             .loc[SAMPLES, "condition"] \
                                             .tolist(),
        gtf=config["reference"]["gtf"]
    run:
        os.makedirs("results/deseq2", exist_ok=True)

        r_script = f"""
        library(DESeq2)
        library(rtracklayer)

        counts <- read.table("{input.counts}",
                             header=TRUE,
                             row.names=1,
                             check.names=FALSE)

        counts <- counts[, -(1:5)]
        colnames(counts) <- sub(".bam$", "", colnames(counts))

        coldata <- data.frame(
            row.names = colnames(counts),
            condition = factor(c({','.join(['"%s"' % c for c in params.conditions])}))
        )

        dds <- DESeqDataSetFromMatrix(
            countData = counts,
            colData   = coldata,
            design    = ~condition
        )

        dds <- dds[rowSums(counts(dds)) > 1, ]
        dds <- DESeq(dds)

        res <- results(dds,
                       contrast=c("condition","L4","D6"))

        # =============================
        # GTFからgene_name取得
        # =============================
        gtf <- import("{params.gtf}")
        gtf_genes <- gtf[gtf$type == "gene"]

        gene_map <- data.frame(
            gene_id = gtf_genes$gene_id,
            gene_name = gtf_genes$gene_name
        )

        res_df <- as.data.frame(res)
        res_df$gene_id <- rownames(res_df)

        res_df <- merge(res_df,
                        gene_map,
                        by="gene_id",
                        all.x=TRUE)

        sig <- subset(res_df,
                      padj < 0.05 &
                      abs(log2FoldChange) > 1)

        write.csv(res_df, file="{output.all_results}", row.names=FALSE)
        write.csv(sig, file="{output.sig_results}", row.names=FALSE)
        saveRDS(dds, file="{output.rld_object}")
        """

        script_path = "results/deseq2/deseq2_tmp.R"
        with open(script_path, "w") as f:
            f.write(r_script)

        shell(f"Rscript {script_path}")
