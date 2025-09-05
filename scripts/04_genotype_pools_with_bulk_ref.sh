#!/bin/bash

#SBATCH --array=1-6
#SBATCH --job-name=pool_genotyping
#SBATCH --account=amc-general
#SBATCH --output=output_pool_genotyping_%A_%a.log
#SBATCH --error=error_pool_genotyping_%A_%a.log
#SBATCH --time=24:00:00
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --mem=2G
#SBATCH --ntasks=10
#SBATCH --nodes=1
#SBATCH --mail-user=grace.akatsu@cuanschutz.edu
#SBATCH --mail-type=ALL 

# Exit immediately on error
set -euo pipefail

##########################################################################################################
# Define paths
##########################################################################################################

ALIGNMENT_DIR="/scratch/alpine/$USER/star_alignment_output"
REF_DIR="/scratch/alpine/$USER/bcftools_genotyping_output_bulk_references"
OUTPUT_DIR="/scratch/alpine/$USER/cellsnp_lite_genotyping_output_pools_BULK_REFS"
mkdir -p "$OUTPUT_DIR"
echo "cellSNP-lite output saving to: $OUTPUT_DIR"

##########################################################################################################
# Gather all BAMs into array
##########################################################################################################

ALL_BAMS=($(find "$ALIGNMENT_DIR" -name "*Aligned.sortedByCoord.out.bam" | sort))
NUM_BAMS=${#ALL_BAMS[@]}
echo "Total BAMs found: $NUM_BAMS"
printf '%s\n' "${ALL_BAMS[@]}"

if [ "$SLURM_ARRAY_TASK_ID" -gt "$NUM_BAMS" ]; then
    echo "Task ID $SLURM_ARRAY_TASK_ID exceeds number of BAMs ($NUM_BAMS). Exiting."
    exit 1
fi

BAM_FILE=${ALL_BAMS[$((SLURM_ARRAY_TASK_ID - 1))]}

# Extract experiment and pool IDs
EXPERIMENT_ID=$(basename "$(dirname "$BAM_FILE")")         
POOL_BAM=$(basename "$BAM_FILE")                                 
POOL_ID=${POOL_BAM%%_*}
EXPERIMENT_POOL_ID="${EXPERIMENT_ID}_${POOL_ID}"

POOL_OUTPUT_DIR="$OUTPUT_DIR/$EXPERIMENT_POOL_ID"
mkdir -p "$POOL_OUTPUT_DIR"

BARCODE_FILE="$ALIGNMENT_DIR/$EXPERIMENT_ID/${POOL_ID}_Solo.out/Gene/filtered/barcodes.tsv"

echo "****** Running cellSNP-lite ******"
echo "Experiment:       $EXPERIMENT_ID"
echo "Pool:             $POOL_ID"
echo "Output:           $POOL_OUTPUT_DIR"
echo "BAM File:         $BAM_FILE"
echo "Barcodes:         $BARCODE_FILE"

##########################################################################################################
# Ensure BAM index exists before running cellSNP-lite
##########################################################################################################

BAM_INDEX="${BAM_FILE}.bai"
if [ ! -f "$BAM_INDEX" ]; then
    echo "Index for $BAM_FILE not found. Creating index..."
    samtools index "$BAM_FILE"
else
    echo "Index for $BAM_FILE already exists."
fi

##########################################################################################################
# Run cellSNP-lite on each pool, using previously generated bulk references for reference
##########################################################################################################

COMPLETION_FILE="$POOL_OUTPUT_DIR/cellsnp_genotype.vcf.gz"

if [ -f "$COMPLETION_FILE" ]; then
    echo "cellSNP-lite output already exists for $EXPERIMENT_POOL_ID. Skipping cellSNP-lite."
    exit 0
fi

cellsnp-lite \
    -s "$BAM_FILE" \
    -O "$POOL_OUTPUT_DIR" \
    -b "$BARCODE_FILE" \
    -R "${REF_DIR}/${EXPERIMENT_POOL_ID}_donor_ref/${EXPERIMENT_POOL_ID}_donor_ref.vcf.gz" \
    -p 20 \
    --minMAF 0.01 \
    --minCOUNT 1 \
    --gzip \
    --cellTAG CB \
    --UMItag UB \
    --genotype

echo "****** Finished cellSNP-lite genotyping on $EXPERIMENT_POOL_ID ******"