#!/bin/bash
#$ -N enrichment
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/5_enrichment/enrichment.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/5_enrichment/enrichment.err
#$ -l h_rt=1:00:00
#$ -l mem_free=5G

# load modules


# load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate basics

# change directories
cd ../../data/gene_enrichment/empiric_counts

# define sets
sets=('altai_2' 'archaic_specific_2' 'chagyrskaya_2' 'denisovan_2' 'introgressed_2' 'neanderthal_2' 'shared_2' 'vindija_2' )

# run enrichment
for s in ${sets[@]}; do python3 ../../../scripts/5_enrichment/enrichment.py --set_name $s; done
echo "complete"
