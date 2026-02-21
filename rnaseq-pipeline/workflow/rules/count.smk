import os
import pandas as pd

# -----------------------------
# サンプル定義
# -----------------------------
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
    run:
        import tempfile

        os.makedirs(os.path.dirname(output.counts), exist_ok=True)

        stranded = config["featurecounts"].get("stranded", 0)
        paired = config["featurecounts"].get("paired", False)

        bam_str = " ".join(input.bams)

        # ---- GTF 解凍対応 ----
        gtf_path = input.gtf
        if gtf_path.endswith(".gz"):
            tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".gtf")
            tmp.close()
            shell(f"gunzip -c {gtf_path} > {tmp.name}")
            gtf_to_use = tmp.name
        else:
            gtf_to_use = gtf_path

        shell(f"""
            featureCounts \
                -T {threads} \
                {'-p' if paired else ''} \
                -s {stranded} \
                -t exon \
                -g gene_id \
                -a {gtf_to_use} \
                -o {output.counts} \
                {bam_str}
        """)

        if gtf_path.endswith(".gz"):
            os.remove(gtf_to_use)


# =============================
# DESeq2
# =============================
rule deseq2:
    input:
        counts="results/counts/gene_counts.txt"
    output:
        all_results="results/deseq2/DESeq2_results_all.csv",
        sig_results="results/deseq2/DESeq2_results_significant.csv",
        rld_object="results/deseq2/rld.rds"
    params:
        # samples.tsv から condition を自動取得（順序も保証）
        conditions=lambda wildcards: samples.set_index("sample")
                                              .loc[SAMPLES, "condition"]
                                              .tolist()
    threads: 1
    run:
        os.makedirs(os.path.dirname(output.all_results), exist_ok=True)

        r_script = f"""
        library(DESeq2)

        counts <- read.table("{input.counts}",
                             header=TRUE,
                             row.names=1,
                             check.names=FALSE)

        # ---- featureCounts のメタ列除去（Chr,Start,End,Strand,Length）----
        counts <- counts[, -(1:5)]

        # ---- 列名から .bam を除去 ----
        colnames(counts) <- sub(".bam$", "", colnames(counts))

        print("=== DEBUG ===")
        print(dim(counts))
        print(colnames(counts))
        print("=============")

        # ---- coldata 作成 ----
        coldata <- data.frame(
            row.names = colnames(counts),
            condition = factor(c({','.join(['"%s"' % c for c in params.conditions])}))
        )

        if(ncol(counts) != nrow(coldata)) {{
            stop("Counts columns and condition length mismatch")
        }}

        # ---- DESeq2 ----
        dds <- DESeqDataSetFromMatrix(
            countData = counts,
            colData   = coldata,
            design    = ~condition
        )

        # 低カウント除去
        dds <- dds[rowSums(counts(dds)) > 1, ]

        dds <- DESeq(dds)

        # L4 vs D6（config.yaml と一致）
        res <- results(dds,
                       contrast=c("condition","L4","D6"))

        sig <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)

        write.csv(as.data.frame(res), file="{output.all_results}")
        write.csv(as.data.frame(sig), file="{output.sig_results}")
        saveRDS(dds, file="{output.rld_object}")
        """

        script_path = os.path.join(
            os.path.dirname(output.all_results),
            "deseq2_tmp.R"
        )

        with open(script_path, "w") as f:
            f.write(r_script)

        shell(f"Rscript {script_path}")

# =============================
# PCA
# =============================
rule pca_plot:
    input:
        rld="results/deseq2/rld.rds"
    output:
        "results/deseq2/DESeq2_PCA.png"
    threads: 1
    run:
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)

        r_script = f"""
        library(DESeq2)
        library(ggplot2)

        dds <- readRDS("{input.rld}")
        rld <- rlog(dds, blind=TRUE)

        pca_data <- plotPCA(rld,
                            intgroup="condition",
                            returnData=TRUE)

        percentVar <- round(100 * attr(pca_data, "percentVar"))

        p <- ggplot(pca_data,
                    aes(PC1, PC2, color=condition)) +
             geom_point(size=3) +
             xlab(paste0("PC1: ", percentVar[1], "% variance")) +
             ylab(paste0("PC2: ", percentVar[2], "% variance")) +
             theme_bw()

        ggsave("{output[0]}", plot=p, width=6, height=5)
        """

        script_path = os.path.join(
            os.path.dirname(output[0]),
            "pca_tmp.R"
        )

        with open(script_path, "w") as f:
            f.write(r_script)

        shell(f"Rscript {script_path}")

