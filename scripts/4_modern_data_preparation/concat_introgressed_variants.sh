#!/bin/bash
#$ -N concat_introgressed_variants
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/4_modern_data_preparation/concat_introgressed_variants.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/4_modern_data_preparation/concat_introgressed_variants.err
#$ -l h_rt=1:00:00
#$ -l mem_free=5G

# load modules

# set path and change directories
cd ../../data/introgression
vernot="../../../../data/vernot_et_al_2016/introgressed_tag_snp_frequencies"

# concatenate Neanderthal haplotypes
cat "$vernot"/all_tag_snps.ASN.merged.ALL.0.3_R2_cluster.1KG_phase3_essentials.bed \
"$vernot"/all_tag_snps.EUR.merged.ALL.0.3_R2_cluster.1KG_phase3_essentials.bed \
"$vernot"/all_tag_snps.PNG.merged.ALL.0.3_R2_cluster.1KG_phase3_essentials.bed \
"$vernot"/all_tag_snps.SAS.merged.ALL.0.3_R2_cluster.1KG_phase3_essentials.bed |
sort -k1,1 -k2,2n > sorted_archaic_introgressed_tag_snps.bed
awk '!seen[$0]++' sorted_archaic_introgressed_tag_snps.bed > sorted_no_dups_archaic_introgressed_tag_snps.bed
