#!/usr/bin/env bash
set -euo pipefail

CONFIG_PATH="config/config.yaml"
CONFIG_DIR=$(dirname "$CONFIG_PATH")
CONFIG_FILE=$(basename "$CONFIG_PATH")
OUT="samples.tsv"

# ヘッダに condition 列を追加
echo -e "sample\tcondition\tlayout\tfq1\tfq2" > "$OUT"

################################
# yq helper
################################
yq_config() {
  docker run --rm \
    -v "$(realpath "$CONFIG_DIR")":/workdir \
    -w /workdir \
    mikefarah/yq "$@" "$CONFIG_FILE"
}

################################
# samples 取得
################################
ALL_SAMPLES=$(yq_config e -r '
  (.samples // {} | keys[]) ,
  (.userdata // {} | keys[])
' 2>/dev/null || true)

for sample in ${ALL_SAMPLES}; do

  # どのブロックか判定
  if yq_config e -r ".samples.${sample}" | grep -vq null; then
    base=".samples.${sample}"
  else
    base=".userdata.${sample}"
  fi

  layout=$(yq_config e -r "${base}.layout")

  # condition 取得
  condition=$(yq_config e -r "${base}.condition")
  if [[ -z "${condition}" || "${condition}" == "null" ]]; then
    echo "[ERROR] condition missing for sample ${sample}"
    exit 1
  fi

  FASTQ_DIR="data/fastq/${sample}"

  ################################
  # SE
  ################################
  if [[ "${layout}" == "SE" ]]; then

    fq1=$(ls ${FASTQ_DIR}/*.fastq.gz 2>/dev/null | head -n 1)

    if [[ -z "${fq1}" ]]; then
      echo "[ERROR] No fastq.gz found for ${sample}"
      exit 1
    fi

    fq2="-"

  ################################
  # PE
  ################################
  elif [[ "${layout}" == "PE" ]]; then

    fq1=$(ls ${FASTQ_DIR}/*_1.fastq.gz 2>/dev/null | head -n 1)
    fq2=$(ls ${FASTQ_DIR}/*_2.fastq.gz 2>/dev/null | head -n 1)

    if [[ -z "${fq1}" || -z "${fq2}" ]]; then
      echo "[ERROR] Paired fastq files missing for ${sample}"
      exit 1
    fi

  else
    echo "[ERROR] layout must be SE or PE"
    exit 1
  fi

  echo -e "${sample}\t${condition}\t${layout}\t${fq1}\t${fq2}" >> "$OUT"

done

echo "[INFO] samples.tsv generated successfully"

################################
# FASTQ を S3 に同期
################################

RNASEQ_S3_BUCKET=$(yq_config e -r '.aws.bucket')
AWS_REGION=$(yq_config e -r '.aws.region')

if [[ -z "$RNASEQ_S3_BUCKET" || "$RNASEQ_S3_BUCKET" == "null" ]]; then
  echo "[ERROR] aws.bucket not set in config"
  exit 1
fi

export RNASEQ_S3_BUCKET AWS_REGION

FASTQ_DIR="data/fastq"

echo "[INFO] Syncing FASTQ directory to S3..."
aws s3 sync "$FASTQ_DIR" "s3://${RNASEQ_S3_BUCKET}/data/fastq/" \
  --exact-timestamps \
  --region "$AWS_REGION"

echo "[INFO] FASTQ S3 sync completed successfully"
