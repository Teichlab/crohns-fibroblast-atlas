#!/bin/bash
#SBATCH -J Distance_Analysis
#SBATCH -p himem
#SBATCH --mem=410G
#SBATCH --time=8:00:00
#SBATCH -o logs/Distances_Xenium.%j.out
#SBATCH -e logs/Distances_Xenium.%j.err

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
module load RStudio-Server/2023.09.1+494-foss-2023a-Java-11-R-4.3.2

#Run Script
echo "[INFO] Running Distance Analysis using Xenium"
Rscript "ULC_Distances_Computation.R"

echo "[INFO] Job end: $(date)"