#!/bin/bash
#$ -N format_fixed_refs
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/other_scripts/format_fixed_refs.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/other_scripts/format_fixed_refs.err
#$ -l h_rt=6:00:00
#$ -l mem_free=30G

# load modules
module load CBI
module load bedtools2/2.30.0
module load htslib/1.13
module load samtools

# change directories
cd ../../data/fixed_refs

# format introgression files
awk '($5=="0") {print $1,$2-1,$2,$3,$4,$5}' OFS='\t' ../introgression/Browning_et_al_2018_Neanderthal_introgressed_variants.txt > formatted_Browning_variants.bed
awk '{print $1,$2,$3,$4,$5,$14}' OFS='\t' ../introgression/Vernot_et_al_2016_introgressed_tag_snps.bed > formatted_Vernot_variants.bed

# intersect
bedtools intersect -a fixed_refs.bed -b formatted_Browning_variants.bed -wa -wb > Browning_variants_intersect.bed
bedtools intersect -a fixed_refs.bed -b formatted_Vernot_variants.bed -wa -wb > Vernot_variants_intersect.bed

# format variants for new FASTA
awk '{print $5,$7,$9}' OFS='\t' Browning_variants_intersect.bed > new_reference/Browning_alt_to_ref_variants.txt
awk '($10!=$9) {print $5,$7,$9}' OFS='\t' Vernot_variants_intersect.bed > new_reference/Vernot_alt_to_ref_variants.txt
cd new_reference
cat Browning_alt_to_ref_variants.txt Vernot_alt_to_ref_variants.txt > tmp.txt
uniq tmp.txt > tmp2.txt
sort -k1,1 -k2n,2 tmp2.txt > alt_to_ref_variants.txt
rm tmp.txt && rm tmp2.txt

# create new FASTA (assuming there are already per chromosome FASTAs present in the working directory)
chrs=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22' )
awk '{print>$1"_variants.txt"}' alt_to_ref_variants.txt
for chr in ${chrs[@]}; do python3 add_variants_to_FASTA.py --variants "$chr"_variants.txt --fasta "$chr".fa --out "$chr"_with_fixed_refs.fa; done
cat chr*_with_fixed_refs.fa > hg19_with_fixed_refs.fa
samtools faidx hg19_with_fixed_refs.fa

# format variants for VCFs
cd ../
awk '{print $5,$7,".",$9,$8,"100","PASS",".","GT","1/1","1/1","1/1","1/1"}' OFS='\t' Browning_variants_intersect.bed > vcfs/Browning_variants.vcf
awk '($10!=$9) {print $5,$7,".",$9,$8,"100","PASS",".","GT","1/1","1/1","1/1","1/1"}' OFS='\t' Vernot_variants_intersect.bed > vcfs/Vernot_variants.vcf
cd vcfs
cat Browning_variants.vcf Vernot_variants.vcf > tmp.vcf
uniq tmp.vcf > tmp2.vcf
sort -k1,1 -k2n,2 tmp2.vcf > combined_unique_variants.vcf
rm tmp.vcf && rm tmp2.vcf

# split files for VCF
split -l 2000 -a 1 combined_unique_variants.vcf fixed_refs_
cat header.txt fixed_refs_a > fixed_refs_1.vcf
cat header.txt fixed_refs_b > fixed_refs_2.vcf
cat header.txt fixed_refs_c > fixed_refs_3.vcf
cat header.txt fixed_refs_d > fixed_refs_4.vcf
cat header.txt fixed_refs_e > fixed_refs_5.vcf
cat header.txt fixed_refs_f > fixed_refs_6.vcf
rm fixed_refs_a && rm fixed_refs_b && rm fixed_refs_c && rm fixed_refs_d && rm fixed_refs_e && rm fixed_refs_f 
