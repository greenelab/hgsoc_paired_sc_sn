#!/bin/sh

#SBATCH --array=1-30
#SBATCH --job-name=star_alignment
#SBATCH --account=amc-general
#SBATCH --output=output_star_alignment_%A_%a.log
#SBATCH --error=error_star_alignment_%A_%a.log
#SBATCH --qos=normal
#SBATCH --time=0-12:00:00
#SBATCH --partition=amilan
#SBATCH --mem=128G
#SBATCH --ntasks=12
#SBATCH --nodes=1
#SBATCH --mail-user=grace.akatsu@cuanschutz.edu
#SBATCH --mail-type=ALL 

set -euo pipefail

# Task-specific ID
export taskID=$SLURM_ARRAY_TASK_ID

# Define paths
PROJ_DIR="${SLURM_SUBMIT_DIR}"
DATA_DIR="/pl/active/cgreene-sc-hgsoc/paired_sc_sn_data_2025"
REF_DATA_DIR="/projects/$USER/hgsoc_paired_sc_sn/reference_data"
MERGED_DIR="/scratch/alpine/$USER/merged_lanes"
mkdir -p "$MERGED_DIR"

OUTPUT_DIR="/scratch/alpine/$USER/star_alignment_output"
mkdir -p "$OUTPUT_DIR"

# Get all unique pool prefixes
ALL_POOLS=($(find "$DATA_DIR"/*R/Fastq -name "*_L*_R1_001.fastq.gz" | xargs -n1 basename | cut -d'_' -f1 | sort | uniq))
POOL_TASK_PREFIX=${ALL_POOLS[$((taskID-1))]}

##########################################################################################################
# Activate conda environment
##########################################################################################################

#echo "Trying to load and activate conda environment..."
#module load miniforge
#source "$(conda info --base)/etc/profile.d/conda.sh"
#conda activate hgsoc_sc_sn_env

##########################################################################################################
# Check if STAR index exists and create if needed
##########################################################################################################

if [ -f "$REF_DATA_DIR/refdata-gex-GRCh38-2024-A/star/SAindex" ]; then
    echo "STAR index exists, skipping."
else
    echo "Creating STAR index..."
    STAR \
        --runMode genomeGenerate \
        --runThreadN 20 \
        --genomeDir "$REF_DATA_DIR/refdata-gex-GRCh38-2024-A/star" \
        --genomeFastaFiles "$REF_DATA_DIR/refdata-gex-GRCh38-2024-A/fasta/genome.fa" \
        --sjdbGTFfile "$REF_DATA_DIR/refdata-gex-GRCh38-2024-A/genes/genes.gtf"
fi

##########################################################################################################
# Run STAR alignment for all matching experiments containing this pool
##########################################################################################################

echo "****** Running STAR alignment on pool: $POOL_TASK_PREFIX ******"

# Identify all experiment directories containing this pool
POOL_FASTQS=($(find "$DATA_DIR"/*R/Fastq -name "${POOL_TASK_PREFIX}_*_L*_R1_001.fastq.gz"))

if [ "${#POOL_FASTQS[@]}" -eq 0 ]; then
    echo "No FASTQ files found for pool $POOL_TASK_PREFIX. Exiting."
    exit 1
fi

for R1_PATH in "${POOL_FASTQS[@]}"; do
    EXPERIMENT_DIR=$(dirname "$R1_PATH")
    EXPERIMENT_ID=$(basename "$(dirname "$EXPERIMENT_DIR")")  # One level up from Fastq dir
    FASTQ_DIR="$EXPERIMENT_DIR"
    OUTPUT_SUBDIR="$OUTPUT_DIR/$EXPERIMENT_ID"
    mkdir -p "$OUTPUT_SUBDIR"

    R1_MERGED="$MERGED_DIR/${EXPERIMENT_ID}_${POOL_TASK_PREFIX}_R1_merged.fastq.gz"
    R2_MERGED="$MERGED_DIR/${EXPERIMENT_ID}_${POOL_TASK_PREFIX}_R2_merged.fastq.gz"

    STAR_OUTPUT_BAM="$OUTPUT_SUBDIR/${POOL_TASK_PREFIX}_Aligned.sortedByCoord.out.bam"

    echo "Processing pool $POOL_TASK_PREFIX in experiment $EXPERIMENT_ID..."

    # Merge lanes
    if [ -f "$R1_MERGED" ] && [ -f "$R2_MERGED" ]; then
        echo "Merged FASTQ files already exist for $POOL_TASK_PREFIX, skipping merge."
    else
        echo "Merging lanes for ${EXPERIMENT_ID}_${POOL_TASK_PREFIX}..."
        # Use the specific FASTQ_DIR identified for this experiment
        cat "$FASTQ_DIR"/${POOL_TASK_PREFIX}_*_L*_R1_001.fastq.gz > "$R1_MERGED"
        cat "$FASTQ_DIR"/${POOL_TASK_PREFIX}_*_L*_R2_001.fastq.gz > "$R2_MERGED"
    fi

    # Run STAR
    if [ -f "$STAR_OUTPUT_BAM" ]; then
        echo "STAR output BAM already exists for $POOL_TASK_PREFIX, skipping alignment."
    else
        if [ $EXPERIMENT_ID = "26384R" ]; then
        # Run alignment on bulk samples
            echo "Running STAR alignment (bulk) for pool $POOL_TASK_PREFIX..."
            cd "$OUTPUT_SUBDIR" || { echo "Failed to cd to $OUTPUT_SUBDIR"; exit 1; }
            STAR \
                --runThreadN 20 \
                --genomeDir "$REF_DATA_DIR/refdata-gex-GRCh38-2024-A/star" \
                --readFilesIn "$R1_MERGED" "$R2_MERGED" \
                --readFilesCommand zcat \
                --outFileNamePrefix "${POOL_TASK_PREFIX}_" \
                --outSAMtype BAM SortedByCoordinate \
                --quantMode TranscriptomeSAM GeneCounts \
                --twopassMode Basic \
                --limitSjdbInsertNsj 3000000
        else
        # Run alignment on sc/sn samples
            echo "Running STAR alignment (sc/sn) for pool $POOL_TASK_PREFIX..."
            
            cd "$OUTPUT_SUBDIR" || { echo "Failed to cd to $OUTPUT_SUBDIR"; exit 1; }
            
            # First pass - run STAR to collect all splice junctions
            STAR \
                --runThreadN 20 \
                --genomeDir "$REF_DATA_DIR/refdata-gex-GRCh38-2024-A/star" \
                --readFilesIn "$R1_MERGED" "$R2_MERGED" \
                --readFilesCommand zcat \
                --outFileNamePrefix "${POOL_TASK_PREFIX}_1stpass_" \
                --outSAMtype None
            
            # Filter SJ.out.tab (only high quality splice junctions)
            echo "Filtering splice junctions..."
            SJ_INPUT="${POOL_TASK_PREFIX}_1stpass_SJ.out.tab"
            SJ_FILTERED="${POOL_TASK_PREFIX}_SJ.filtered.tab"
            
            # Only keep junctions supported by >= 5 uniquely mapped reads
            # and with maximum overhang of >=20 base pairs
            awk '$7 >= 5 && $9 >= 20' "$SJ_INPUT" > "$SJ_FILTERED"
            
            # Second pass - run STAR with filtered junctions
            STAR \
                --runThreadN 20 \
                --genomeDir "$REF_DATA_DIR/refdata-gex-GRCh38-2024-A/star" \
                --readFilesIn "$R1_MERGED" "$R2_MERGED" \
                --readFilesCommand zcat \
                --sjdbFileChrStartEnd "$SJ_FILTERED" \
                --outFileNamePrefix "${POOL_TASK_PREFIX}_" \
                --outSAMtype BAM SortedByCoordinate \
                --quantMode TranscriptomeSAM GeneCounts \
                --limitBAMsortRAM 100000000000 \
                --soloType CB_UMI_Simple \
                --soloCBstart 1 --soloCBlen 16 \
                --soloUMIstart 17 --soloUMIlen 12 \
                --soloBarcodeReadLength 0 \
                --soloCBwhitelist None \
                --outSAMattributes NH HI AS nM CB UB
        fi
    fi

    echo "Finished pool $POOL_TASK_PREFIX in experiment $EXPERIMENT_ID."
done

echo "****** All alignments completed for $POOL_TASK_PREFIX ******"