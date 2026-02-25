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
# helper: yq via docker
################################
yq_config() {
    docker run --rm \
        -v "$CONFIG_DIR":/workdir \
        -w /workdir \
        mikefarah/yq "$@" "$CONFIG_FILE"
}

################################
# 1. AWS config
################################
RNASEQ_S3_BUCKET=$(yq_config e -r '.aws.bucket')
AWS_REGION=$(yq_config e -r '.aws.region')

if [[ -z "$RNASEQ_S3_BUCKET" || "$RNASEQ_S3_BUCKET" == "null" ]]; then
    echo "ERROR: aws.bucket is not set in config"
    exit 1
fi

export RNASEQ_S3_BUCKET AWS_REGION

################################
# 2. system packages (Ubuntu 24.04対応)
################################
echo "[INFO] Installing system packages..."
sudo apt update
sudo apt install -y \
    wget \
    curl \
    unzip \
    build-essential \
    libncurses-dev \
    zlib1g-dev \
    pigz

################################
# 3. AWS CLI v2 install (公式方式)
################################
echo "[INFO] Installing AWS CLI v2..."
cd /tmp

curl -s "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip -q awscliv2.zip
sudo ./aws/install --update

rm -rf aws awscliv2.zip

echo "[INFO] AWS CLI version:"
aws --version

################################
# SRA Toolkit install
################################
sudo apt install -y sra-toolkit
################################
# 動作確認
################################
echo "[INFO] Checking SRA Toolkit..."
fasterq-dump --version

################################
# 5. check tools
################################
for tool in fasterq-dump pigz aws; do
    command -v $tool >/dev/null 2>&1 || { echo "[ERROR] $tool not found"; exit 1; }
done
echo "[INFO] Tools installed successfully"

################################
# 6. S3 bucket / prefixes 初期化
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

################################
# 7. Create working directories
################################
RAW_DIR="$WORKDIR/data/fastq"
REF_FASTA_DIR="$WORKDIR/data/reference/fasta"
REF_GTF_DIR="$WORKDIR/data/reference/gtf"

mkdir -p "$RAW_DIR" "$RAW_DIR/tmp" "$REF_FASTA_DIR" "$REF_GTF_DIR"

################################
# 8. Download Trimmomatic adapters
################################
ADAPTERS_DIR="$WORKDIR/adapters"
mkdir -p "$ADAPTERS_DIR"

echo "[INFO] Downloading Trimmomatic adapters..."

SE_ADAPTER_URL="https://github.com/usadellab/Trimmomatic/raw/master/adapters/TruSeq3-SE.fa"
SE_ADAPTER_FILE="$ADAPTERS_DIR/TruSeq3-SE.fa"

if [ ! -f "$SE_ADAPTER_FILE" ]; then
    wget -q "$SE_ADAPTER_URL" -O "$SE_ADAPTER_FILE"
    echo "[INFO] SE adapter downloaded"
fi

PE_ADAPTER_URL="https://github.com/usadellab/Trimmomatic/raw/master/adapters/TruSeq3-PE.fa"
PE_ADAPTER_FILE="$ADAPTERS_DIR/TruSeq3-PE.fa"

if [ ! -f "$PE_ADAPTER_FILE" ]; then
    wget -q "$PE_ADAPTER_URL" -O "$PE_ADAPTER_FILE"
    echo "[INFO] PE adapter downloaded"
fi

echo
echo "[INFO] Environment setup completed!"

