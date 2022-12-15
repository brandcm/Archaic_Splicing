#!/bin/bash
#$ -N format_introgressed_refs
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/other_scripts/format_introgressed_refs.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/other_scripts/format_introgressed_refs.err
#$ -l h_rt=6:00:00
#$ -l mem_free=30G

# load modules
module load CBI
module load bedtools2/2.30.0
module load htslib/1.13
module load samtools

# change directories
cd ../../data/introgressed_refs

# format variants for new FASTA
#awk '{print $1,$2,$4}' OFS='\t' introgressed_ref_tag.txt > new_reference/alt_to_ref_variants.txt

# create new FASTA (assuming there are already per chromosome FASTAs present in the working directory)
cd new_reference
chrs=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22' )
awk '{print>$1"_variants.txt"}' alt_to_ref_variants.txt
for chr in ${chrs[@]}; do python3 add_variants_to_FASTA.py --variants "$chr"_variants.txt --fasta "$chr".fa --out "$chr"_with_introgressed_refs.fa; done
cat chr*_with_introgressed_refs.fa > hg19_with_introgressed_refs.fa
samtools faidx hg19_with_introgressed_refs.fa

# format variants for VCFs
cd ../
awk '{print $1,$2,".",$4,$3,"100","PASS",".","GT","1/1","1/1","1/1","1/1"}' OFS='\t' introgressed_ref_tag.txt > vcfs/introgressed_refs.vcf

# split files for VCF
# I used a larger split here, recommend doing a smaller one and adding the appropriate code below
cd vcfs
split -l 5000 -a 1 introgressed_refs.vcf introgressed_refs_
cat header.txt introgressed_refs_a > introgressed_refs_1.vcf
cat header.txt introgressed_refs_b > introgressed_refs_2.vcf
rm introgressed_refs_a && rm introgressed_refs_b
