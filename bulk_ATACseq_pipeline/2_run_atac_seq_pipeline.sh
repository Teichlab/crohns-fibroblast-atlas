#!/bin/bash

# to run:
# conda activate atac_env
# bsub < 2_run_atac_seq_pipeline.sh
# 
# installation instructions: https://github.com/ENCODE-DCC/atac-seq-pipeline

# LSF configuration
#BSUB -J ENCODE_ATAC_seq[1-3]  # should match number of .json files in config directory
#BSUB -n 6                #number of cores
#BSUB -q long             # queue
#BSUB -R "select[mem>32G] rusage[mem=32G] span[hosts=1]"
#BSUB -M 32G
#BSUB -o logs/output_%J_%I.out #output. %J is job-id %I is job-array index
#BSUB -e logs/error_%J_%I.err

mkdir -p logs

# Load software
module load /software/modules/cellgen/singularity


# nth ls file in config directory
# full path to .json config file
input_json=$(ls -d $PWD/config/*.json | sed -n "$LSB_JOBINDEX"p)

# Print path of .json config file used
echo $input_json

# Create output directory, changing working directory such that output files are written here by the subsequent caper run
# basename without .json extension
mkdir -p results/$(basename -s .json $input_json)
cd results/$(basename -s .json $input_json)

# Run ENCODE ATAC-seq pipeline locally with caper
caper init local
caper run /nfs/team205/sk29/software/atac-seq-pipeline/atac.wdl -i "${input_json}" --singularity  #--max-concurrent-tasks 1