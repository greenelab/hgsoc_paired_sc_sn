#!/bin/sh

#SBATCH --array=1-6
#SBATCH --job-name=genotype_bulk_reference_in_pools
#SBATCH --account=amc-general
#SBATCH --output=output_bulk_ref_in_pools_%A_%a.log
#SBATCH --error=error_bulk_ref_in_pools_%A_%a.log
#SBATCH --time=09:00:00
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --mem=16G
#SBATCH --ntasks=10
#SBATCH --nodes=1
#SBATCH --mail-user=grace.akatsu@cuanschutz.edu
#SBATCH --mail-type=ALL 

set -euo pipefail

##########################################################################################################
# Paths
##########################################################################################################
PROJ_DIR="${SLURM_SUBMIT_DIR}"
ALIGNMENT_DIR="/scratch/alpine/$USER/star_alignment_output/26384R"
OUTPUT_DIR="/scratch/alpine/$USER/cellsnp_lite_genotyping_output_bulk_references_per_sample"
mkdir -p "$OUTPUT_DIR"
REF_SNPS="/projects/$USER/hgsoc_paired_sc_sn/reference_data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz"

##########################################################################################################
# Define pools and bulk BAMs
##########################################################################################################
POOLS=("26210R_26210X1" "26210R_26210X2" "26224R_26224X1" "26224R_26224X2" "26225R_26225X1" "26225R_26225X2")

BULK_BAMS=(
    "26384X15_Aligned.sortedByCoord.out.bam 26384X14_Aligned.sortedByCoord.out.bam 26384X16_Aligned.sortedByCoord.out.bam 26384X20_Aligned.sortedByCoord.out.bam"
    "26384X17_Aligned.sortedByCoord.out.bam 26384X13_Aligned.sortedByCoord.out.bam 26384X19_Aligned.sortedByCoord.out.bam 26384X23_Aligned.sortedByCoord.out.bam"
    "26384X18_Aligned.sortedByCoord.out.bam 26384X21_Aligned.sortedByCoord.out.bam 26384X24_Aligned.sortedByCoord.out.bam"
    "26384X19_Aligned.sortedByCoord.out.bam 26384X23_Aligned.sortedByCoord.out.bam 26384X18_Aligned.sortedByCoord.out.bam 26384X21_Aligned.sortedByCoord.out.bam"
    "26384X15_Aligned.sortedByCoord.out.bam 26384X16_Aligned.sortedByCoord.out.bam 26384X17_Aligned.sortedByCoord.out.bam 26384X13_Aligned.sortedByCoord.out.bam"
    "26384X18_Aligned.sortedByCoord.out.bam 26384X21_Aligned.sortedByCoord.out.bam 26384X24_Aligned.sortedByCoord.out.bam"
)

##########################################################################################################
# Select pool and loop through individual BAMs
##########################################################################################################
POOL_ID="${POOLS[$((SLURM_ARRAY_TASK_ID-1))]}"
BAMS="${BULK_BAMS[$((SLURM_ARRAY_TASK_ID-1))]}"
POOL_OUT="$OUTPUT_DIR/${POOL_ID}_donor_ref"
mkdir -p "$POOL_OUT"

VCF_LIST=()
for bam in $BAMS; do
    FULL_BAM="$ALIGNMENT_DIR/$bam"
    SAMPLE_ID="$(basename "$bam" .bam | cut -d'_' -f1)"
    SAMPLE_OUT="$OUTPUT_DIR/individual_bulk_genotypes/${SAMPLE_ID}_genotype"
    mkdir -p "$SAMPLE_OUT"

    SAMPLE_VCF="$SAMPLE_OUT/cellSNP.cells.vcf.gz"

    # Check if sample genotyping already done
    if [ -f "$SAMPLE_VCF" ]; then
        echo "Sample $SAMPLE_ID already genotyped, skipping..."
    else
        # Check BAM index
        if [ ! -f "${FULL_BAM}.bai" ]; then
            echo "BAM index not found for $FULL_BAM, indexing now..."
            samtools index "$FULL_BAM"
        fi

        echo "Running cellSNP-lite for $FULL_BAM -> $SAMPLE_OUT"

        cellsnp-lite \
            -s "$FULL_BAM" \
            -O "$SAMPLE_OUT" \
            -R "$REF_SNPS" \
            -p 10 \
            --cellTAG None \
            --UMItag None \
            --minMAF 0.05 \
            --minCOUNT 1 \
            --gzip \
            --genotype
    fi

    # Add sample VCF to merge list
    VCF_LIST+=("$SAMPLE_VCF")
done


# Merge all sample VCFs into one donor reference VCF for the pool
POOL_VCF="$POOL_OUT/${POOL_ID}_donor_ref.vcf.gz"

if [ -f "$POOL_VCF" ]; then
    echo "Pool-level donor reference VCF for $POOL_ID already exists, skipping merge."
else
    echo "Merging individual sample VCFs for pool $POOL_ID ..."
    bcftools merge -m all -Oz -o "$POOL_VCF" "${VCF_LIST[@]}"
    bcftools index "$POOL_VCF"
    echo "Finished donor reference VCF for pool $POOL_ID -> $POOL_VCF"
fi

