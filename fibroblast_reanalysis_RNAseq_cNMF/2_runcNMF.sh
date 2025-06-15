#!/bin/bash
# Submits cNMF jobs in parallel on farm.
# runtime: ~20m
#
# conda activate cnmf_env

# number of workers, split for each k to evaluate
n_workers=20

# Raw counts h5ad file
h5ad_file="../data/rna/scRNA_batch1_2.raw.fibroblasts.h5ad"

# Normalize gene counts and prepare run parameters
# -c path to h5ad file
# -k space separated list of K values to be tested
# --numgenes high variance genes to be used

cnmf prepare \
    --output-dir ./results \
    --name fibroblast_cNMF \
    -c "$h5ad_file" \
    --n-iter 100 \
    -k {5..25} \
    --seed 14 \
    --total-workers "$n_workers" \
    --numgenes 2000

mkdir -p logs  # folder for LSF logs, if not already

# Loop over jobs, running cnmf factorize
for (( i = 0; i < $n_workers; i++ )); do
    bsub -G teichlab -q normal \
        -n 2 \
        -M 2GB -R "select[mem>2GB] rusage[mem=2GB]" \
        -o logs/output.log -e logs/error.log \
        cnmf factorize \
            --output-dir ./results \
            --name fibroblast_cNMF \
            --worker-index "$i" \
            --total-workers "$n_workers"
done