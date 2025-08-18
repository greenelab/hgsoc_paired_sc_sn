#!/bin/bash

#SBATCH --array=1-6
#SBATCH --job-name=vireo_demultiplex
#SBATCH --account=amc-general
#SBATCH --output=output_vireo_demultiplex_%A_%a.log
#SBATCH --error=error_vireo_demultiplex_%A_%a.log
#SBATCH --time=18:00:00
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --mem=24G
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mail-user=grace.akatsu@cuanschutz.edu
#SBATCH --mail-type=ALL 

# Exit immediately on error
set -euo pipefail

################################################################################
# Define paths
################################################################################

# Set base directory where all pools are located
BASE_DIR="/scratch/alpine/$USER/cellsnp_lite_genotyping_output_pools"

# Reference directory where the genotypes of the reference bulks in pools are located
REF_DIR="/scratch/alpine/$USER/cellsnp_lite_genotyping_output_bulk_references"

# Set output dir
VIREO_OUT_DIR="/scratch/alpine/$USER/vireo_genetic_demultiplexing_output"
    mkdir -p "$VIREO_OUT_DIR"

################################################################################
# Set up
################################################################################

# Number of donors in each pool
N_DONOR=4

echo "Starting Vireo demultiplexing for all pools in: $BASE_DIR"
echo "Expecting $N_DONOR donors per pool"
echo "---------------------------------------------"

################################################################################
# Get list of pool directories and select one based on SLURM_ARRAY_TASK_ID
################################################################################

# Create an array of all pool subdirectories
mapfile -t POOL_DIRS < <(find "$BASE_DIR" -mindepth 1 -maxdepth 1 -type d | sort)

# Check if index is valid
if (( SLURM_ARRAY_TASK_ID > ${#POOL_DIRS[@]} )); then
    echo "Error: SLURM_ARRAY_TASK_ID ($SLURM_ARRAY_TASK_ID) exceeds number of pools (${#POOL_DIRS[@]})"
    exit 1
fi

# Select the pool directory for this task
POOL_DIR="${POOL_DIRS[$((SLURM_ARRAY_TASK_ID-1))]}"
POOL_NAME=$(basename "$POOL_DIR")
VIREO_OUT_DIR_POOL="${VIREO_OUT_DIR}/${POOL_NAME}"
mkdir -p "$VIREO_OUT_DIR_POOL"

################################################################################
# Run Vireo
################################################################################

echo "Processing $POOL_NAME..."
echo "Expecting $N_DONOR donors per pool"

# Sanity check: make sure cellSNP output exists
if [ ! -f "${POOL_DIR}/cellSNP.tag.AD.mtx" ]; then
    echo "  Skipping $POOL_NAME — cellSNP output not found."
    exit 0
fi

vireo \
  --cellData="$POOL_DIR" \
  --donorFile="${REF_DIR}/${POOL_NAME}_donor_ref.vcf.gz" \
  --nDonor="$N_DONOR" \
  --outDir="$VIREO_OUT_DIR_POOL"

echo "✅ Done with $POOL_NAME → results saved to $VIREO_OUT_DIR_POOL"
