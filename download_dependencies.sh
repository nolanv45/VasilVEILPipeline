#!/bin/bash
set -e

BASE_URL="https://zenodo.org/records/18113859/files"

mkdir -p tools/embeddings containers/umap containers/phidra

download() {
    local file=$1
    local out=$2

    if [ ! -f "$out" ]; then
        echo "Downloading $file"
        wget -c "$BASE_URL/$file" -O "$out"
    else
        echo "$out already exists"
    fi
}

# ESM Model
download "esm2_t36_3B_UR50D.pt" tools/embeddings/esm2_t36_3B_UR50D.pt
download "3b_model_checkpoint.pt" tools/embeddings/3b_model_checkpoint.pt

# Singularity containers
download "umap.sif" containers/umap/umap.sif
download "phidra.sif" containers/phidra/phidra.sif

echo "All dependencies downloaded from Zenodo."