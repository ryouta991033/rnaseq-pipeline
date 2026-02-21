#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("Usage: run_deseq2.R <counts.tsv> <metadata.tsv> <config.yaml> <output.tsv>")
}

count_file  <- args[1]
meta_file   <- args[2]
config_file <- args[3]
out_file    <- args[4]

suppressPackageStartupMessages({
  library(DESeq2)
  library(yaml)
})

############################
# load config
############################
config <- yaml::read_yaml(config_file)

if (is.null(config$deseq2)) {
  stop("deseq2 section is missing in config.yaml")
}

deseq_cfg <- config$deseq2

# required keys
required_keys <- c("design", "contrast")
missing <- setdiff(required_keys, names(deseq_cfg))
if (length(missing) > 0) {
  stop(paste("Missing deseq2 config keys:", paste(missing, collapse = ", ")))
}

if (is.null(deseq_cfg$contrast$factor) ||
    is.null(deseq_cfg$contrast$numerator) ||
    is.null(deseq_cfg$contrast$denominator)) {
  stop("deseq2.contrast must contain factor, numerator, denominator")
}

design_str <- deseq_cfg$design
factor     <- deseq_cfg$contrast$factor
num        <- deseq_cfg$contrast$numerator
den        <- deseq_cfg$contrast$denominator

############################
# load counts (featureCounts)
############################
counts_raw <- read.delim(
  count_file,
  comment.char = "#",
  check.names = FALSE
)

if (!"Geneid" %in% colnames(counts_raw)) {
  stop("counts file must contain 'Geneid' column (featureCounts output)")
}

rownames(counts_raw) <- counts_raw$Geneid
counts <- counts_raw[, -(1:6)]

############################
# load metadata
############################
coldata <- read.delim(meta_file, check.names = FALSE)

if (!"sample" %in% colnames(coldata)) {
  stop("metadata.tsv must contain a column named 'sample'")
}

rownames(coldata) <- coldata$sample

############################
# sanity checks
############################
if (!all(colnames(counts) %in% rownames(coldata))) {
  stop("Sample names in counts do not match metadata")
}

counts <- counts[, rownames(coldata)]

if (!factor %in% colnames(coldata)) {
  stop(paste("Factor not found in metadata:", factor))
}

############################
# reference level (optional)
############################
if (!is.null(deseq_cfg$reference)) {
  for (factor_name in names(deseq_cfg$reference)) {
    ref_level <- deseq_cfg$reference[[factor_name]]

    if (!factor_name %in% colnames(coldata)) {
      stop(paste("Reference factor not found in metadata:", factor_name))
    }

    coldata[[factor_name]] <- relevel(
      as.factor(coldata[[factor_name]]),
      ref = ref_level
    )
  }
}

############################
# DESeq2
############################
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = coldata,
  design    = as.formula(design_str)
)

dds <- DESeq(dds)

############################
# contrast
############################
res <- results(
  dds,
  contrast = c(factor, num, den)
)

############################
# output
############################
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

write.table(
  res_df,
  file = out_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

