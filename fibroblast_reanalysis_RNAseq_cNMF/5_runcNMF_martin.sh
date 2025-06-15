#!/bin/bash
# Submits cNMF jobs in parallel on farm.
# runtime: ~20m
#
# conda activate cnmf_env

# number of workers, split for each k to evaluate
n_workers=24

# Raw counts h5ad file
# h5ad_file="/nfs/team205/re5/data/martin19_all.processed.cellxgene_for_rasa.h5ad"
h5ad_file="/lustre/scratch126/cellgen/team205/sk29/matthias_fb/data/external/martin19_raw/martin19_all.raw.annot.h5ad"

# Normalize gene counts and prepare run parameters
# -c path to h5ad file
# -k space separated list of K values to be tested
# --numgenes high variance genes to be used

cnmf prepare \
    --output-dir ./results \
    --name martin_cNMF \
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
        -o logs/martin_output.log -e logs/martin_error.log \
        cnmf factorize \
            --output-dir ./results \
            --name martin_cNMF \
            --worker-index "$i" \
            --total-workers "$n_workers"
done