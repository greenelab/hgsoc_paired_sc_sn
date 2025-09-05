#!/bin/sh

#SBATCH --job-name=fastq_file_preprocessing
#SBATCH --account=amc-general
#SBATCH --output=output_fastq_file_preprocessing.log
#SBATCH --error=error_fastq_file_preprocessing.log
#SBATCH --time=12:00:00
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --ntasks-per-node=32
#SBATCH --nodes=1 

# Exit if any command fails
set -e

# Define paths
PROJ_DIR="${SLURM_SUBMIT_DIR}"
DATA_DIR="/pl/active/cgreene-sc-hgsoc/paired_sc_sn_data_2025"
OUTPUT_DIR="/projects/$USER/hgsoc_paired_sc_sn/preprocessing_fastp_qc_info"

# Ensure OUTPUT_DIR exists
if [ ! -d "$OUTPUT_DIR" ]; then
    mkdir -p "$OUTPUT_DIR"
    echo "Created directory: $OUTPUT_DIR"
fi

##########################################################################################################
# Activate conda environment
##########################################################################################################

module load miniforge

# Check if 'hgsoc_sc_sn_env' exists; if not, create it and include fastp
if ! conda env list | grep -q "hgsoc_sc_sn_env"; then
    echo "Creating conda environment 'hgsoc_sc_sn_env'..."
    conda create -y -n hgsoc_sc_sn_env -c bioconda -c conda-forge fastp star cellsnp-lite vireo
fi

# Activate environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate hgsoc_sc_sn_env

##########################################################################################################
# Run fastp quality control on raw reads data (FASTQ files)
##########################################################################################################

echo "****** Running fastp quality control on raw reads data ******"

find "$DATA_DIR" -type f -path "*/Fastq/*_R1_001.fastq.gz" | while read -r R1; do
    SAMPLE=$(basename "$R1" _R1_001.fastq.gz)
    R2=$(echo "$R1" | sed 's/_R1_001.fastq.gz/_R2_001.fastq.gz/')
    
    REPORT_HTML="$OUTPUT_DIR/${SAMPLE}_fastp_report.html"
    REPORT_JSON="$OUTPUT_DIR/${SAMPLE}_fastp_report.json"

    if [ -f "$REPORT_HTML" ] && [ -f "$REPORT_JSON" ]; then
        echo "Skipping sample $SAMPLE: already processed."
        continue
    fi

    echo "Processing sample: $SAMPLE"

    TMP_OUT1=$(mktemp)
    TMP_OUT2=$(mktemp)

    fastp \
        -i "$R1" -I "$R2" \
        -o "$TMP_OUT1" -O "$TMP_OUT2" \
        --html "$REPORT_HTML" \
        --json "$REPORT_JSON" \
        --thread 8

done

echo "****** fastp processing complete. Results saved in: $OUTPUT_DIR ******"
