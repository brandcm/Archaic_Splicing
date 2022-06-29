#!/bin/bash
#$ -N empiric_FDR
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/5_enrichment/empiric_FDR.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/5_enrichment/empiric_FDR.err
#$ -l h_rt=30:00:00
#$ -l mem_free=10G

# load modules


# load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate basics

# change directories
cd ../../data/gene_enrichment/empiric_counts

# define sets
sets=('altai_2' 'archaic_specific_2' 'chagyrskaya_2' 'denisovan_2' 'introgressed_2' 'neanderthal_2' 'shared_2' 'vindija_2' )

# run FDR script
for s in ${sets[@]}; do for o in GWAS HPO; do python3 ../../../scripts/5_enrichment/empiric_FDR.py --input ${s}_${o}_empiric_counts.txt --subset 1000 --output ../empiric_FDR/${s}_${o}_empiric_FDR.txt; done; done
echo "complete"
