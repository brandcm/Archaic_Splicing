#!/bin/bash
#$ -N get_1KG_sample_names
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/analysis/get_1KG_sample_names.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/analysis/get_1KG_sample_names.err
#$ -l h_rt=1:00:00
#$ -l mem_free=5G

# load modules
module load CBI
module load bcftools/1.13

# change directories
cd ../../data/thousand_genomes

# run script
bcftools query -l subset_concat_normalized_1000_genomes.vcf.gz > 1KG_sample_names.txt
