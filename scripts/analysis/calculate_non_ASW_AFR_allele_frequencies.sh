#!/bin/bash
#$ -N calculate_non_ASW_AFR_allele_frequencies
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/5_analysis/calculate_non_ASW_AFR_allele_frequencies.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/5_analysis/calculate_non_ASW_AFR_allele_frequencies.err
#$ -l h_rt=6:00:00
#$ -l mem_free=20G

# load modules
module load CBI
module load bcftools/1.13

# change directories
cd ../../data/archaic_variants_in_humans

# get any archaic variants also present in 1000 Genomes, splice altering or not
bcftools view  -S ../thousand_genomes/non_ASW_AFR_1000_genomes_samples.txt archaic_variants_subset_concat_normalized_1000_genomes.vcf.gz -O z -o archaic_variants_subset_concat_normalized_1000_genomes_non_ASW_AFR.vcf.gz
bcftools index -t archaic_variants_subset_concat_normalized_1000_genomes_non_ASW_AFR.vcf.gz

# query chrom, pos, ref allele, alt allele, allele count, allele number, and allele frequencies
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\n' archaic_variants_subset_concat_normalized_1000_genomes_non_ASW_AFR.vcf.gz > non_ASW_AFR_allele_frequencies_hg38.txt
echo "query complete"
