#!/bin/bash
#$ -N all_archaic_sites
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/other_scripts/all_archaic_sites.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/other_scripts/all_archaic_sites.err
#$ -l h_rt=60:00:00
#$ -l mem_free=60G

# load modules
module load CBI
module load bcftools/1.13

# set path to hg19 and change directories
hg19="../../../../data/hg19_fasta/2022-03-14/hg19.fa"
cd ../../data/archaic_genomes

# list chromosomes
chrs=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22' 'chrX' )

# normalize VCFs but split combined records into biallelic records and filter out sites with same ref allele and low quality sites & genotypes
for c in ${chrs[@]}; do bcftools norm -m - -f "$hg19" combined_archaic_"$c".vcf.gz | bcftools filter -i 'QUAL>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools view -O z -o all_filtered_combined_archaic_"$c".vcf.gz & done
wait
for c in ${chrs[@]}; do bcftools index -t all_filtered_combined_archaic_"$c".vcf.gz & done
wait

# concat and sort
bcftools concat all_filtered_combined_archaic_chr*.vcf.gz -O z -o all_filtered_combined_archaic_snvs_2.vcf.gz
bcftools sort all_filtered_combined_archaic_snvs_2.vcf.gz -O z -o all_filtered_combined_archaic_snvs.vcf.gz
bcftools index -t all_filtered_combined_archaic_snvs.vcf.gz
rm all_filtered_combined_archaic_snvs_2.vcf.gz

# get genic sites and query
for c in  ${chrs[@]}; do bcftools view -R ../annotations/grch37_gene_annotations.bed all_filtered_combined_archaic_"$c".vcf.gz | bcftools view -O z -o genic_all_"$c".vcf.gz; done
for c in  ${chrs[@]}; do bcftools index -t genic_all_"$c".vcf.gz; done
bcftools concat genic_all_chr*.vcf.gz -O z -o concat_genic_all.vcf.gz
bcftools index -t concat_genic_all.vcf.gz

rm genic_all_chr*.vcf.gz

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' concat_genic_all.vcf.gz > ../fixed_refs/concat_genic_all.txt
awk '{print $1,$2-1,$2,$3}' OFS='\t' ../fixed_refs/concat_genic_all.txt > ../fixed_refs/fixed_refs.bed
