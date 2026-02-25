#!/usr/bin/env bash
set -euo pipefail

################################
# Usage
################################
if [[ $# -lt 1 ]]; then
  echo "Usage: $0 [restore|backup] [--delete]"
  exit 1
fi

MODE="$1"
DELETE_FLAG=""
[[ "${2:-}" == "--delete" ]] && DELETE_FLAG="--delete"

################################
# 固定パス
################################
WORKDIR=$(realpath "./")
CONFIG_PATH="$WORKDIR/config/config.yaml"
CONFIG_DIR=$(dirname "$CONFIG_PATH")
CONFIG_FILE=$(basename "$CONFIG_PATH")

################################
# yq helper（docker経由）
################################
yq_config() {
  docker run --rm \
    -v "$CONFIG_DIR":/workdir \
    -w /workdir \
    mikefarah/yq "$@" "$CONFIG_FILE"
}

################################
# config読み込み
################################
BUCKET=$(yq_config e -r '.aws.bucket')
REGION=$(yq_config e -r '.aws.region')
PREFIX=$(yq_config e -r '.aws.results_prefix')

[[ -z "$BUCKET" || "$BUCKET" == "null" ]] && { echo "aws.bucket not set"; exit 1; }
[[ -z "$PREFIX" || "$PREFIX" == "null" ]] && { echo "aws.prefix not set"; exit 1; }

export AWS_DEFAULT_REGION="$REGION"

################################
# Local paths
################################
LOCAL_DATA="${WORKDIR}/data"
LOCAL_RESULTS="${WORKDIR}/results"

mkdir -p "$LOCAL_DATA" "$LOCAL_RESULTS"

################################
# S3 paths（現在の構造に合わせる）
################################
# dataはbucket直下
S3_DATA="s3://${BUCKET}/data"

# resultsはprefix配下
S3_RESULTS="s3://${BUCKET}/${PREFIX}/results"

################################
# Sync logic
################################
if [[ "$MODE" == "restore" ]]; then

  echo "======================================"
  echo "[INFO] RESTORE MODE (S3 → EBS)"
  echo "======================================"

  echo "[INFO] Restoring DATA from $S3_DATA"
  aws s3 sync "$S3_DATA" "$LOCAL_DATA" $DELETE_FLAG

  echo "[INFO] Restoring RESULTS from $S3_RESULTS"
  aws s3 sync "$S3_RESULTS" "$LOCAL_RESULTS" $DELETE_FLAG

elif [[ "$MODE" == "backup" ]]; then

  echo "======================================"
  echo "[INFO] BACKUP MODE (EBS → S3)"
  echo "======================================"

  echo "[INFO] Backing up DATA to $S3_DATA"
  aws s3 sync "$LOCAL_DATA" "$S3_DATA" $DELETE_FLAG

  echo "[INFO] Backing up RESULTS to $S3_RESULTS"
  aws s3 sync "$LOCAL_RESULTS" "$S3_RESULTS" $DELETE_FLAG

else
  echo "ERROR: mode must be restore or backup"
  exit 1
fi

echo
echo "======================================"
echo "[INFO] SYNC COMPLETED SUCCESSFULLY"
echo "======================================"

