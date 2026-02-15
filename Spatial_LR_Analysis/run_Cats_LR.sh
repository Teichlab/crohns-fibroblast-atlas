#!/bin/bash
#SBATCH -J LR_Analysis
#SBATCH -p short
#SBATCH --cpus-per-task=6
#SBATCH --mem=80G
#SBATCH --time=4:00:00
#SBATCH -o logs/Cats_LR.%j.out  # stdout
#SBATCH -e logs/Cats_LR.%j.err  # stderr

#Information of the job
echo "------------------------------------------------"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "------------------------------------------------"

set -euo pipefail

#Load software modules required
echo "[INFO] Loading required modules..."
module purge
module load R/4.5.1-gfbf-2023a-bare-noSciPy

#Run Merge Script
echo "[INFO] Running L-R Analysis using Xenium"
Rscript Niche_9_11_LR_Analysis_CatsCraddle.R

echo "[INFO] Job end: $(date)"