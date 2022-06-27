#!/bin/bash
#$ -N concat_spliceai_snv_vcfs
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/concat_spliceai_snv_vcfs.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/concat_spliceai_snv_vcfs.err
#$ -l h_rt=6:00:00
#$ -l mem_free=20G

# load modules
module load CBI
module load bcftools/1.13

# concatenate files
#bcftools concat ../data/split.x*.vcf -O z -o ../data/splice_annotated_archaic_snvs_2.vcf.gz  
#bcftools sort ../data/splice_annotated_archaic_snvs_2.vcf.gz -O z -o ../data/splice_annotated_archaic_snvs.vcf.gz
#echo "concat complete"
#rm ../data/splice_annotated_archaic_snvs_2.vcf.gz
#rm ../data/split.x*.vcf
#rm ../data/split.x*.vcf.gz 
#rm ../data/split.x*.vcf.gz.tbi

# index concantenated file
#bcftools index -t ../data/splice_annotated_archaic_snvs.vcf.gz
#echo "indexing complete"

# get relevant fields, remove everything in INFO that's not the alternate allele, turn pipes into tabs, turn commas into new lines, and fill in CHROM, POS, and REF for extra annotations
bcftools query -f '%CHROM\t%POS\t%REF\t[%GT\t]%INFO\n' ../data/splice_annotated_archaic_snvs.vcf.gz | awk '{sub(/.*SpliceAI=/, "", $8)}1' | tr '|' '\t' | tr ',' '\n' |
awk 'NF==10{print p"\t"q"\t"r"\t"s"\t"t"\t"u"\t"v"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10; next} {p=$1}{q=$2}{r=$3}{s=$4}{t=$5}{u=$6}{v=$7} 1' |
awk 'NF==17{print}{}' | tr ' ' '\t' > ../data/spliceai_archaic_snvs.txt

# delete error file if empty
if [ ! -s concat_spliceai_snv_vcfs.err ] ; then
  rm concat_spliceai_snv_vcfs.err
fi