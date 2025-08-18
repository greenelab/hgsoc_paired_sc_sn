#!/bin/sh

#SBATCH --array=1-6
#SBATCH --job-name=pool_genotyping
#SBATCH --account=amc-general
#SBATCH --output=output_pool_genotyping_%A_%a.log
#SBATCH --error=error_pool_genotyping_%A_%a.log
#SBATCH --time=12:00:00
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --mem=5G
#SBATCH --ntasks=10
#SBATCH --nodes=1
#SBATCH --mail-user=grace.akatsu@cuanschutz.edu
#SBATCH --mail-type=ALL 

# Exit immediately on error
set -euo pipefail

##########################################################################################################
# Define paths
##########################################################################################################

PROJ_DIR="${SLURM_SUBMIT_DIR}"
ALIGNMENT_DIR="/scratch/alpine/$USER/star_alignment_output"
OUTPUT_DIR="/scratch/alpine/$USER/cellsnp_lite_genotyping_output_pools"
mkdir -p "$OUTPUT_DIR"
echo "cellSNP-lite output saving to: $OUTPUT_DIR"

##########################################################################################################
# Ensure SNP reference has "chr" prefixes to match BAMs
##########################################################################################################

FIXED_REF_SNPS="/projects/$USER/hgsoc_paired_sc_sn/reference_data/genome1K.phase3.SNP_AF5e2.hg38.chr.vcf.gz"

if [ ! -f "$FIXED_REF_SNPS" ]; then
    echo "Creating hg38-compatible SNP reference with chr prefixes..."
    bcftools annotate \
        --rename-chrs <(awk '{print $1"\tchr"$1}' <(seq 1 22; echo X; echo Y; echo MT)) \
        -Oz -o "$FIXED_REF_SNPS" \
        "/projects/$USER/hgsoc_paired_sc_sn/reference_data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz" \
    && tabix -p vcf "$FIXED_REF_SNPS" \
    || { echo "Error: failed to create $FIXED_REF_SNPS"; exit 1; }
else
    echo "Found existing hg38-compatible SNP reference: $FIXED_REF_SNPS"
fi

REF_SNPS="$FIXED_REF_SNPS"

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
echo "Reference SNPs:   $REF_SNPS"

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
# Run cellSNP-lite on each pool
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
    -R "$REF_SNPS" \
    -p 20 \
    --minMAF 0.05 \
    --minCOUNT 1 \
    --gzip \
    --cellTAG CB \
    --UMItag UB \
    --genotype

echo "****** Finished cellSNP-lite genotyping on $EXPERIMENT_POOL_ID ******"