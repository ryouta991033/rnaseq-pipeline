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

AWS_BUCKET=$(yq_config e -r '.aws.bucket')
AWS_REGION=$(yq_config e -r '.aws.region')
AWS_PREFIX=$(yq_config e -r '.aws.results_prefix // ""')

if [[ -z "$AWS_BUCKET" || "$AWS_BUCKET" == "null" ]]; then
  echo "[ERROR] aws.bucket not set in config"
  exit 1
fi

BASE_PATH="s3://${AWS_BUCKET}"
if [[ -n "$AWS_PREFIX" && "$AWS_PREFIX" != "null" ]]; then
  BASE_PATH="${BASE_PATH}/${AWS_PREFIX}"
fi

echo "[INFO] Checking S3 results path:"
echo "       ${BASE_PATH}"
echo

################################
# results ディレクトリ
################################
echo "====================================="
echo " Results directory"
echo "====================================="
aws s3 ls "${BASE_PATH}/results/" \
  --recursive \
  --human-readable \
  --summarize \
  --region "$AWS_REGION" || echo "[WARN] results not found"
echo

################################
# snakemake_logs
################################
echo "====================================="
echo " Snakemake logs"
echo "====================================="
aws s3 ls "${BASE_PATH}/snakemake_logs/" \
  --recursive \
  --human-readable \
  --summarize \
  --region "$AWS_REGION" || echo "[WARN] snakemake_logs not found"
echo

################################
# run log
################################
echo "====================================="
echo " snakemake_run.log"
echo "====================================="
aws s3 ls "${BASE_PATH}/snakemake_run.log" \
  --region "$AWS_REGION" || echo "[WARN] snakemake_run.log not found"
echo

################################
# stats file
################################
echo "====================================="
echo " snakemake_stats.json"
echo "====================================="
aws s3 ls "${BASE_PATH}/snakemake_stats.json" \
  --region "$AWS_REGION" || echo "[WARN] snakemake_stats.json not found"
echo

echo "[INFO] S3 results check completed"

