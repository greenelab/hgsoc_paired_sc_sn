# Analysis of paired single cell RNA-seq and single nucleus RNA-seq data from HGSOC human samples
## Scripts:
- *01_preprocess_data.sh*
    This script uses fastp to run QC on sequencing quality.
- *02_star_alignment.sh*
    This script uses STAR and STARsolo to align the single cell, single nucleus, and bulk RNA-seq data to a reference genome.
- *03_genotype_pools.sh*
    This script uses cellSNP-lite to genotype the pooled single cell and single nucleus RNA-seq data.
- *04_genotype_bulk_references_in_pools.sh*
    This script uses cellSNP-lite to genotype the bulk RNA-seq data and combines this genotyping into a pool to use as reference for genetic demultiplexing.
- *05_vireo_genetic_demultiplex.sh*
    This script uses vireo to genetically demultiplex the pooled single cell and single nucleus RNA-seq data.

## Reference data:
### Reference SNP list:
https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz

### Reference genome:
https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz
