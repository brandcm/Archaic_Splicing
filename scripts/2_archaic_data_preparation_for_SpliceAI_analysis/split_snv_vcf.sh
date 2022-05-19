#!/bin/bash
#$ -N split_snv_vcf
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/2_archaic_data_preparation_for_spliceai_analysis/split_snv_vcf.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/2_archaic_data_preparation_for_spliceai_analysis/split_snv_vcf.err
#$ -l h_rt=12:00:00
#$ -l mem_free=30G

# load modules
module load CBI
module load bcftools/1.13

# change directories
cd ../../data/archaic_genomes

# split combined VCF into smaller files by 5000 SNVs
bcftools query -f'%CHROM\t%POS\n' concat_subset_filtered_combined_archaic_snvs.vcf.gz | split -l 5000 -a 3 --numeric-suffixes=100
for file in x*; do bcftools view -T $file concat_subset_filtered_combined_archaic_snvs.vcf.gz -O z -o split_$file.vcf.gz & done
wait
echo "split complete"

# index
for f in split_x*.vcf.gz; do bcftools index -t $f & done
wait
echo "indexing complete"

# remove SNP files used to generate split VCFs
# MAKE SURE NO OTHER FILES BEGIN WITH "x"
rm x*
