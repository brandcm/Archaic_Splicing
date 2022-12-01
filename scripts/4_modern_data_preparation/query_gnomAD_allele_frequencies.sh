#!/bin/bash
#$ -N query_gnomAD_allele_frequencies
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/4_modern_data_preparation/query_gnomAD_allele_frequencies.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/4_modern_data_preparation/query_gnomAD_allele_frequencies.err
#$ -l h_rt=24:00:00
#$ -l mem_free=60G

# load modules
module load CBI
module load bcftools/1.13

# change directories
cd ../../data/archaic_variants_in_humans/

# set path
gnomAD_vcfs=../../../../data/wynton_databases/gnomAD/3.0/variants/genomes

# set variable
chrs=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22' 'chrX' )

# get any archaic variants also present in gnomAD, splice altering or not
for chr in ${chrs[@]}; do bcftools view -f 'PASS,.' -R ../archaic_variants/all_variant_sites_hg38.bed "$gnomAD_vcfs"/gnomad.genomes.r3.0.sites."$chr".vcf.bgz -O z -o "$chr"_archaic_variants_in_gnomAD.vcf.gz; done
for chr in ${chrs[@]}; do bcftools index -t "$chr"_archaic_variants_in_gnomAD.vcf.gz; done

# concat
bcftools concat *_archaic_variants_in_gnomAD.vcf.gz -O z -o archaic_variants_in_gnomAD.vcf.gz
bcftools index -t archaic_variants_in_gnomAD.vcf.gz

# query chrom, pos, ref allele, alt allele, allele count, allele number, and variant_type
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\t%AF\t%variant_type\n' archaic_variants_in_gnomAD.vcf.gz > archaic_variants_in_gnomAD_allele_frequencies.txt
echo "query complete"
