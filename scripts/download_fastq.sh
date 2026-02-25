#!/usr/bin/env bash
set -euo pipefail

################################
# 固定パス設定
################################
WORKDIR=$(realpath "./")
CONFIG_PATH="$WORKDIR/config/config.yaml"
CONFIG_DIR=$(dirname "$CONFIG_PATH")
CONFIG_FILE=$(basename "$CONFIG_PATH")

echo "[INFO] Using config: $CONFIG_PATH"
echo "[INFO] Working directory: $WORKDIR"
echo

################################
# helper: yq
################################
yq_config() {
    docker run --rm \
        -v "$CONFIG_DIR":/workdir \
        -w /workdir \
        mikefarah/yq "$@" "$CONFIG_FILE"
}

################################
# AWS config
################################
RNASEQ_S3_BUCKET=$(yq_config e -r '.aws.bucket')
AWS_REGION=$(yq_config e -r '.aws.region')
if [[ -z "$RNASEQ_S3_BUCKET" || "$RNASEQ_S3_BUCKET" == "null" ]]; then
    echo "ERROR: aws.bucket is not set in config"
    exit 1
fi
export RNASEQ_S3_BUCKET AWS_REGION

################################
# S3 初期化
################################
echo "[INFO] Initializing S3 bucket and prefixes..."
if ! aws s3 ls "s3://$RNASEQ_S3_BUCKET" >/dev/null 2>&1; then
    echo "[INFO] Bucket $RNASEQ_S3_BUCKET does not exist. Creating..."
    aws s3 mb "s3://$RNASEQ_S3_BUCKET" --region "$AWS_REGION"
fi

PREFIXES=("data/fastq/" "data/reference/fasta/" "data/reference/gtf/")
for prefix in "${PREFIXES[@]}"; do
    if ! aws s3 ls "s3://$RNASEQ_S3_BUCKET/$prefix" >/dev/null 2>&1; then
        echo "[INFO] Creating S3 prefix: $prefix"
        aws s3api put-object --bucket "$RNASEQ_S3_BUCKET" --key "$prefix"
    fi
done
echo "[INFO] S3 initialization done"
echo

################################
# paths
################################
RAW_DIR="$WORKDIR/data/fastq"
REF_FASTA_DIR="$WORKDIR/data/reference/fasta"
REF_GTF_DIR="$WORKDIR/data/reference/gtf"
mkdir -p "$RAW_DIR" "$RAW_DIR/tmp" "$REF_FASTA_DIR" "$REF_GTF_DIR"

################################
# threads
################################
THREADS_FASTERQ=$(yq_config e -r '.threads.fasterq')
THREADS_PIGZ=$(yq_config e -r '.threads.pigz')

################################
# layout consistency check
################################
LAYOUTS=$(yq_config e -r '.samples[].layout' | sort -u)
if [[ $(echo "$LAYOUTS" | wc -l) -ne 1 ]]; then
    echo "ERROR: mixed single/paired samples are not supported"
    exit 1
fi
PROJECT_LAYOUT=$(echo "$LAYOUTS" | head -1)
echo "[INFO] project layout = $PROJECT_LAYOUT"

################################
# reference download
################################
FASTA_URL=$(yq_config e -r '.download.reference.fasta.url')
FASTA_NAME=$(yq_config e -r '.download.reference.fasta.name')
GTF_URL=$(yq_config e -r '.download.reference.gtf.url')
GTF_NAME=$(yq_config e -r '.download.reference.gtf.name')

download_with_curl() {
    local url="$1"
    local out="$2"
    [[ -f "$out" ]] && echo "[INFO] skip existing $(basename "$out")" && return
    echo "[INFO] downloading $(basename "$out")"
    curl -L "$url" -o "$out"
}

download_with_curl "$FASTA_URL" "$REF_FASTA_DIR/$FASTA_NAME"
download_with_curl "$GTF_URL" "$REF_GTF_DIR/$GTF_NAME"

################################
# FASTQ download (並列化対応)
################################
SAMPLES=$(yq_config e -r '.samples | keys | .[]')
TMP_LIST=$(mktemp)

for sample in $SAMPLES; do
    SRA=$(yq_config e -r ".samples.${sample}.sra")
    LAYOUT=$(yq_config e -r ".samples.${sample}.layout")
    echo "$sample $SRA $LAYOUT" >> "$TMP_LIST"
done

CPU_PARALLEL=2  # 同時サンプル数

cat "$TMP_LIST" | xargs -P "$CPU_PARALLEL" -I {} bash -c '
set -- {}
SAMPLE=$1
SRA=$2
LAYOUT=$3
RAW_DIR="'"$RAW_DIR"'"
THREADS_FASTERQ="'"$THREADS_FASTERQ"'"
THREADS_PIGZ="'"$THREADS_PIGZ"'"

SAMPLE_DIR="$RAW_DIR/$SAMPLE"
mkdir -p "$SAMPLE_DIR"
mkdir -p "$RAW_DIR/tmp"

echo "[INFO] downloading $SAMPLE ($SRA) [$LAYOUT]"

if [[ "$LAYOUT" == "paired" ]]; then
    fasterq-dump "$SRA" --split-files -O "$SAMPLE_DIR" -e "$THREADS_FASTERQ" --temp "$RAW_DIR/tmp"
    pigz -p "$THREADS_PIGZ" "$SAMPLE_DIR"/*.fastq
    mv "$SAMPLE_DIR/${SRA}_1.fastq.gz" "$SAMPLE_DIR/${SAMPLE}_R1.fastq.gz"
    mv "$SAMPLE_DIR/${SRA}_2.fastq.gz" "$SAMPLE_DIR/${SAMPLE}_R2.fastq.gz"
else
    fasterq-dump "$SRA" -O "$SAMPLE_DIR" -e "$THREADS_FASTERQ" --temp "$RAW_DIR/tmp"
    pigz -p "$THREADS_PIGZ" "$SAMPLE_DIR"/*.fastq
    mv "$SAMPLE_DIR/${SRA}.fastq.gz" "$SAMPLE_DIR/${SAMPLE}.fastq.gz"
fi

echo "[INFO] uploading $SAMPLE to S3"
aws s3 sync "$SAMPLE_DIR" "s3://'"$RNASEQ_S3_BUCKET"'/data/fastq/$SAMPLE"

echo "[INFO] done $SAMPLE"
' _

rm -f "$TMP_LIST"

################################
# reference upload (ループ外でまとめて)
################################
echo "[INFO] uploading reference FASTA directory to S3"
aws s3 sync "$REF_FASTA_DIR" "s3://$RNASEQ_S3_BUCKET/data/reference/fasta/" --exact-timestamps

echo "[INFO] uploading reference GTF directory to S3"
aws s3 sync "$REF_GTF_DIR" "s3://$RNASEQ_S3_BUCKET/data/reference/gtf/" --exact-timestamps


echo "[INFO] all FASTQ downloads & uploads completed successfully"

