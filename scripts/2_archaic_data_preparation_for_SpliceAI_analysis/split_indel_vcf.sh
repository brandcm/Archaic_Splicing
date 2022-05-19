#!/bin/bash
#$ -N split_indel_vcf
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/2_archaic_data_preparation_for_spliceai_analysis/split_indel_vcf.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/2_archaic_data_preparation_for_spliceai_analysis/split_indel_vcf.err
#$ -l h_rt=12:00:00
#$ -l mem_free=30G

# load modules
module load CBI
module load bcftools/1.13

# change directories
cd ../../data/archaic_genomes

# split combined VCF into smaller files by 5000 InDels
bcftools query -f'%CHROM\t%POS\n' subset_filtered_combined_archaic_indels.vcf.gz > indel_splits.txt
split -l 5000 -a 3 --numeric-suffixes=100 indel_splits.txt a
for file in a*; do bcftools view -T $file subset_filtered_combined_archaic_indels.vcf.gz -O z -o split_indel_$file.vcf.gz & done
wait
echo "split indels complete"

# index
for f in split_indel_a*.vcf.gz; do bcftools index -t $f & done
wait
echo "indexing indels complete"

# remove files used to generate split VCFs
# MAKE SURE NO OTHER FILES BEGIN WITH "a"
rm indel_splits.txt
rm a*
