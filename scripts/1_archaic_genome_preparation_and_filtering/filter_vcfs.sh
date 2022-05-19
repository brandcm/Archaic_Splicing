#!/bin/bash
#$ -N filter_vcfs
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/1_archaic_genome_preparation_and_filtering/filter_vcfs.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/1_archaic_genome_preparation_and_filtering/filter_vcfs.err
#$ -l h_rt=12:00:00
#$ -l mem_free=30G

# load modules
module load CBI
module load bcftools/1.13

# set path to hg19 and change directories
hg19="../../../../data/hg19_fasta/2022-03-14/hg19.fa"
cd ../../data/archaic_genomes

# list chromosomes
chrs=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22' 'chrX' )

# normalize VCFs but split combined records into biallelic records and filter out sites with same ref allele and low quality sites & genotypes
for c in ${chrs[@]}; do bcftools norm -m - -f "$hg19" combined_archaic_"$c".vcf.gz | bcftools filter -e 'ALT="."' | bcftools filter -i 'QUAL>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools view -O z -o filtered_combined_archaic_"$c".vcf.gz & done
wait
for c in ${chrs[@]}; do bcftools index -t filtered_combined_archaic_"$c".vcf.gz & done
wait

bcftools norm -m - -f "hg19" combined_archaic_indels.vcf.gz | bcftools filter -e 'ALT="."' | bcftools filter -i 'QUAL>=30' | bcftools filter -S . -i 'FMT/DP>=10 & FMT/GQ>=30' | bcftools view -O z -o filtered_combined_archaic_indels.vcf.gz
bcftools index -t filtered_combined_archaic_indels.vcf.gz

echo "filtering and indexing complete"
