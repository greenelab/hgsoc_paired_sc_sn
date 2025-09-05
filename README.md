
# Analysis of paired single cell, single nucleus, and bulk RNA-seq data from HGSOC human samples
## Scripts:
### **01_preprocess_data.sh**
This script uses fastp to run QC on sequencing quality.
### **02_star_alignment.sh**
This script uses STAR and STARsolo to align the single cell, single nucleus, and bulk RNA-seq data to a reference genome.
### **03_genotype_pools.sh**
This script uses cellSNP-lite to genotype the pooled single cell and single nucleus RNA-seq data. It uses a set of reference SNPs from the alignment reference genome.
### **04_genotype_bulk_references_in_pools.sh**
This script uses cellSNP-lite to genotype the bulk RNA-seq data and combines this genotyping into a pool to use as reference for genetic demultiplexing.
### **05_vireo_genetic_demultiplex.sh**
This script uses vireo to genetically demultiplex the pooled single cell and single nucleus RNA-seq data.

## Reference data:
### Reference genome:
https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz
Unzip the reference genome files and place in directory reference_data, on the same level as directory scripts

