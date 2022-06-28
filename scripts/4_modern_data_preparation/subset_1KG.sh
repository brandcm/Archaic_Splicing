#!/bin/bash
#$ -N subset_1KG
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/4_modern_data_preparation/subset_1KG.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/4_modern_data_preparation/subset_1KG.err
#$ -l h_rt=12:00:00
#$ -l mem_free=30G

# load modules
module load CBI
module load bcftools/1.13

# designate path to 1000 genomes data, hg19 fasta, and change directories
thousand_genomes=../../../../data/wynton_databases/1000_genomes/release/20190312_biallelic_SNV_and_INDEL/
hg38_fasta=../../../../data/hg38_fasta/2022-03-14/hg38.fa
cd ../../data/thousand_genomes

# list chromosomes
chrs=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22' 'chrX' )

# normalize 1000 genomes  
for c in ${chrs[@]}; do bcftools norm -m - -f "$hg38_fasta" "$thousand_genomes"ALL."$c".shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz -O z -o normalized_"$c"_1000_genomes.vcf.gz & done
wait
echo "normalization complete"

# index files
for c in ${chrs[@]}; do bcftools index -t normalized_"$c"_1000_genomes.vcf.gz & done
wait
echo "indexing complete"

# concat, sort, and index
bcftools concat normalized_chr*_1000_genomes.vcf.gz -O z -o concat_normalized_1000_genomes.vcf.gz
bcftools sort concat_normalized_1000_genomes_2.vcf.gz -O z -o concat_normalized_1000_genomes.vcf.gz
bcftools index -t concat_normalized_1000_genomes.vcf.gz

# subset to genes 
bcftools view -R ../annotations/grch38_gene_annotations.bed concat_normalized_1000_genomes.vcf.gz -O z -o subset_concat_normalized_1000_genomes.vcf.gz
bcftools index -t subset_concat_normalized_1000_genomes.vcf.gz
echo "subset complete"

# remove files that are no longer needed
rm normalized_chr*_1000_genomes.vcf.gz
rm normalized_chr*_1000_genomes.vcf.gz.tbi
rm concat_normalized_1000_genomes.vcf.gz
rm concat_normalized_1000_genomes.vcf.gz.tbi
echo "temp files removed"
