#!/bin/bash
#$ -N concat_introgressed_refs
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/other_scripts/concat_introgressed_refs.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/other_scripts/concat_introgressed_refs.err
#$ -l h_rt=6:00:00
#$ -l mem_free=20G

# load modules
module load CBI
module load bcftools/1.13

cd ../../data/introgressed_refs/vcfs

# concat files
# multiple tail commands results in a header so let's remove
# also given the append, let's check to make sure the file doesn't exist
if [ ! -f introgressed_ref_annotations.txt ]; then
    tail -n +29 introgressed_refs_*_annotated.vcf | grep '^[^=]' >> introgressed_ref_annotations.txt
fi

# get relevant fields, remove everything in INFO that's not the alternate allele, turn pipes into tabs, turn commas into new lines, and fill in CHROM, POS, and REF for extra annotations
awk '{print $1,$2,$4,$10,$11,$12,$13,$8}' OFS='\t' introgressed_ref_annotations.txt | awk '{sub(/.*SpliceAI=/, "", $8)}1' | tr '|' '\t' | tr ',' '\n' |
awk 'NF==10{print p"\t"q"\t"r"\t"s"\t"t"\t"u"\t"v"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10; next} {p=$1}{q=$2}{r=$3}{s=$4}{t=$5}{u=$6}{v=$7} 1' |
awk 'NF==17{print}{}' | tr ' ' '\t' > ../spliceai_annotated_introgressed_refs.txt
