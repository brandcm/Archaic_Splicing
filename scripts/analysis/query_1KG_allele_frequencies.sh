#!/bin/bash
#$ -N query_1KG_allele_frequencies
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/analysis/query_1KG_allele_frequencies.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/analysis/query_1KG_allele_frequencies.err
#$ -l h_rt=12:00:00
#$ -l mem_free=30G

# load modules
module load CBI
module load bcftools/1.13

# change directories
cd ../../data/archaic_variants_in_humans/

# get any archaic variants also present in 1000 Genomes, splice altering or not
bcftools view -R ../archaic_variants/all_variant_sites_hg38_no_chr.bed ../thousand_genomes/subset_concat_normalized_1000_genomes.vcf.gz -O z -o archaic_variants_subset_concat_normalized_1000_genomes.vcf.gz
bcftools index -t archaic_variants_subset_concat_normalized_1000_genomes.vcf.gz

# query chrom, pos, ref allele, alt allele, allele count, allele number, and allele frequencies
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\t%AF\t%EAS_AF\t%EUR_AF\t%AFR_AF\t%AMR_AF\t%SAS_AF\n' archaic_variants_subset_concat_normalized_1000_genomes.vcf.gz > allele_frequencies_hg38.txt
echo "query complete"
