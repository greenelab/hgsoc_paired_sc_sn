#!/bin/sh

#SBATCH --array=1-6
#SBATCH --job-name=genotype_bulk_reference_in_pools
#SBATCH --account=amc-general
#SBATCH --output=output_bulk_ref_in_pools_genotyping_%A_%a.log
#SBATCH --error=error_bulk_ref_in_pools_genotyping_%A_%a.log
#SBATCH --time=09:00:00
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --mem=16G
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
ALIGNMENT_DIR="/scratch/alpine/$USER/star_alignment_output/26384R" # Experiment 26384R has the bulks
OUTPUT_DIR="/scratch/alpine/$USER/cellsnp_lite_genotyping_output_bulk_references_in_pools"
mkdir -p "$OUTPUT_DIR"
REF_SNPS="/projects/$USER/hgsoc_paired_sc_sn/reference_data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz"

echo "cellSNP-lite output saving genotyping of bulk references in pools to: $OUTPUT_DIR"

##########################################################################################################
# Associate bulk samples with pools
##########################################################################################################

# Map SLURM_ARRAY_TASK_ID to pool ID and list of BAMs for that pool
# Mappings can be found in /pl/active/cgreene-sc-hgsoc/paired_sc_sn_data_2025/scSeqOvaryR01-20250703_summary_file.xlsx
# There are corresponding dissociated bulk and tissue bulk available to use as reference
# Here I will uses tissue bulks, as one patient sample does not have a dissociated (cell) bulk
# Note that one sample that has neither tissue nor dissociated bulk reference!

POOLS=("26210R_26210X1" "26210R_26210X2" "26224R_26224X1" "26224R_26224X2" "26225R_26225X1" "26225R_26225X2")

# Each string in BULK_BAMS contains space-separated BAM paths for that pool
BULK_BAMS=(
    "26384X15_Aligned.sortedByCoord.out.bam 26384X14_Aligned.sortedByCoord.out.bam 26384X16_Aligned.sortedByCoord.out.bam 26384X20_Aligned.sortedByCoord.out.bam"   # 26210R_26210X1
    "26384X17_Aligned.sortedByCoord.out.bam 26384X13_Aligned.sortedByCoord.out.bam 26384X19_Aligned.sortedByCoord.out.bam 26384X23_Aligned.sortedByCoord.out.bam"   # 26210R_26210X2 - missing one bulk reference!
    "26384X18_Aligned.sortedByCoord.out.bam 26384X21_Aligned.sortedByCoord.out.bam 26384X24_Aligned.sortedByCoord.out.bam"                                          # 26224R_26224X1
    "26384X19_Aligned.sortedByCoord.out.bam 26384X23_Aligned.sortedByCoord.out.bam 26384X18_Aligned.sortedByCoord.out.bam 26384X21_Aligned.sortedByCoord.out.bam"   # 26224R_26224X2
    "26384X15_Aligned.sortedByCoord.out.bam 26384X16_Aligned.sortedByCoord.out.bam 26384X17_Aligned.sortedByCoord.out.bam 26384X13_Aligned.sortedByCoord.out.bam"   # 26225R_26225X1
    "26384X18_Aligned.sortedByCoord.out.bam 26384X21_Aligned.sortedByCoord.out.bam 26384X24_Aligned.sortedByCoord.out.bam"                                          # 26225R_26225X2  - missing one bulk reference!
)

##########################################################################################################
# Select pool and merge BAMs
##########################################################################################################

POOL_ID="${POOLS[$((SLURM_ARRAY_TASK_ID-1))]}"
BAMS="${BULK_BAMS[$((SLURM_ARRAY_TASK_ID-1))]}"

# Prepend ALIGNMENT_DIR to each BAM filename
MERGED_BAM_LIST=""
for bam in $BAMS; do
    FULL_BAM="$ALIGNMENT_DIR/$bam"
    if [ ! -f "$FULL_BAM" ]; then
        echo "Warning: BAM $FULL_BAM not found, skipping"
        continue
    fi
    MERGED_BAM_LIST="$MERGED_BAM_LIST $FULL_BAM"
done

if [ -z "$MERGED_BAM_LIST" ]; then
    echo "Error: No BAMs found for pool $POOL_ID. Exiting."
    exit 1
fi

MERGED_BAM="$OUTPUT_DIR/${POOL_ID}_bulk_merged.bam"

if [ -f "$MERGED_BAM" ] && [ -f "${MERGED_BAM}.bai" ]; then
    echo "Merged BAM already exists for $POOL_ID. Skipping merge."
else
    echo "Merging bulk BAMs for pool $POOL_ID: $MERGED_BAM_LIST"
    samtools merge -@ 10 -o "$MERGED_BAM" $MERGED_BAM_LIST
    samtools index "$MERGED_BAM"
fi

##########################################################################################################
# Run cellSNP-lite genotyping
##########################################################################################################

POOL_OUT="$OUTPUT_DIR/${POOL_ID}_pooled_bulks_genotype"
mkdir -p "$POOL_OUT"

cellsnp-lite \
    -s "$MERGED_BAM" \
    -O "$POOL_OUT" \
    -R "$REF_SNPS" \
    -p 10 \
    --cellTAG None \
    --UMItag None \
    --minMAF 0.05 \
    --minCOUNT 1 \
    --gzip \
    --genotype \

echo "Finished cellSNP-lite genotyping on $POOL_OUT"

