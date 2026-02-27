################################
# plots.smk
################################

import os

# =============================
# Volcano plot（Top10遺伝子名表示）
# =============================
rule volcano_plot:
    input:
        rld="results/deseq2/rld.rds",
        res="results/deseq2/DESeq2_results_all.csv"
    output:
        "results/plots/volcano_plot.pdf"
    run:
        os.makedirs("results/plots", exist_ok=True)

        r_script = f"""
        library(ggplot2)
        library(ggrepel)

        res_df <- read.csv("{input.res}")

        res_df$significant <- "NS"
        res_df$significant[
            res_df$padj < 0.05 &
            abs(res_df$log2FoldChange) > 1
        ] <- "Significant"

        topgenes <- head(res_df[order(res_df$padj), ], 10)

        p <- ggplot(res_df,
                    aes(x=log2FoldChange,
                        y=-log10(padj),
                        color=significant)) +
             geom_point(alpha=0.6) +
             geom_text_repel(data=topgenes,
                             aes(label=gene_name),
                             size=3) +
             theme_bw()

        ggsave("{output[0]}", plot=p, width=6, height=5)
        """

        script_path = "results/plots/volcano_tmp.R"
        with open(script_path, "w") as f:
            f.write(r_script)

        shell(f"Rscript {script_path}")


# =============================
# Heatmap（gene_name表示）
# =============================
rule heatmap:
    input:
        rld="results/deseq2/rld.rds",
        res="results/deseq2/DESeq2_results_all.csv"
    output:
        "results/plots/heatmap.pdf"
    run:
        os.makedirs("results/plots", exist_ok=True)

        r_script = f"""
        library(DESeq2)
        library(pheatmap)

        dds <- readRDS("{input.rld}")
        rld <- rlog(dds, blind=TRUE)

        res_df <- read.csv("{input.res}")
        res_df <- res_df[!is.na(res_df$padj), ]

        sig <- subset(res_df,
                      padj < 0.05 &
                      abs(log2FoldChange) > 1)

        sig <- sig[order(sig$log2FoldChange,
                         decreasing=TRUE), ]

        topgenes <- sig$gene_id[1:min(50, nrow(sig))]

        mat <- assay(rld)[topgenes, , drop=FALSE]
        rownames(mat) <- sig$gene_name[match(topgenes,
                                             sig$gene_id)]

        mat <- mat - rowMeans(mat)

        pdf("{output[0]}", width=6, height=8)
        pheatmap(mat,
                 cluster_rows=FALSE,
                 cluster_cols=TRUE,
                 show_rownames=TRUE)
        dev.off()
        """

        script_path = "results/plots/heatmap_tmp.R"
        with open(script_path, "w") as f:
            f.write(r_script)

        shell(f"Rscript {script_path}")
