#!/usr/bin/env bash
set -euo pipefail

################################
# config 読み込み
################################
CONFIG_PATH="config/config.yaml"
CONFIG_DIR=$(dirname "$CONFIG_PATH")
CONFIG_FILE=$(basename "$CONFIG_PATH")

yq_config() {
  docker run --rm \
    -v "$(realpath "$CONFIG_DIR")":/workdir \
    -w /workdir \
    mikefarah/yq "$@" "$CONFIG_FILE"
}

RNASEQ_S3_BUCKET=$(yq_config e -r '.aws.bucket')
AWS_REGION=$(yq_config e -r '.aws.region')

if [[ -z "$RNASEQ_S3_BUCKET" || "$RNASEQ_S3_BUCKET" == "null" ]]; then
  echo "[ERROR] aws.bucket not set in config"
  exit 1
fi

echo "[INFO] Checking S3 bucket: s3://${RNASEQ_S3_BUCKET}"
echo

################################
# FASTQ
################################
echo "=============================="
echo " FASTQ files (data/fastq)"
echo "=============================="
aws s3 ls "s3://${RNASEQ_S3_BUCKET}/data/fastq/" \
  --recursive \
  --human-readable \
  --summarize \
  --region "$AWS_REGION"

echo

################################
# Reference FASTA
################################
echo "=============================="
echo " Reference FASTA"
echo "=============================="
aws s3 ls "s3://${RNASEQ_S3_BUCKET}/data/reference/fasta/" \
  --recursive \
  --human-readable \
  --summarize \
  --region "$AWS_REGION"

echo

################################
# Reference GTF
################################
echo "=============================="
echo " Reference GTF"
echo "=============================="
aws s3 ls "s3://${RNASEQ_S3_BUCKET}/data/reference/gtf/" \
  --recursive \
  --human-readable \
  --summarize \
  --region "$AWS_REGION"

echo
echo "[INFO] S3 check completed"

