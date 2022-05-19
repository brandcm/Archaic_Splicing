#!/bin/bash
#$ -N subset_vcfs
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/2_archaic_data_preparation_for_spliceai_analysis/subset_vcfs.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/2_archaic_data_preparation_for_spliceai_analysis/subset_vcfs.err
#$ -l h_rt=12:00:00
#$ -l mem_free=30G

# load modules
module load CBI
module load bcftools/1.13

# change directories
cd ../../data/archaic_genomes

# list chromosomes
chrs=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22' 'chrX' )

# subset VCFs to InDels and SNVs within genes
for c in ${chrs[@]}; do bcftools view -R ../annotations/slopped_grch37_gene_annotations.bed filtered_combined_archaic_$c.vcf.gz | bcftools view -O z -o subset_filtered_combined_archaic_$c.vcf.gz & done
wait

bcftools view -R ../annotations/slopped_grch37_gene_annotations.bed filtered_combined_archaic_indels.vcf.gz | bcftools view -O z -o subset_filtered_combined_archaic_indels.vcf.gz
echo "subset complete"

# index subset
for c in ${chrs[@]}; do bcftools index -t subset_filtered_combined_archaic_$c.vcf.gz & done
wait

bcftools index -t subset_filtered_combined_archaic_indels.vcf.gz
echo "subset indexing complete"

# concatenate SNV files
bcftools concat_subset_filtered_combined_archaic_chr*.vcf.gz -O z -o concat_subset_filtered_combined_archaic_snvs_2.vcf.gz  
bcftools sort concat_subset_filtered_combined_archaic_snvs_2.vcf.gz -O z -o concat_subset_filtered_combined_archaic_snvs.vcf.gz
echo "merge and sort complete"
rm concat_subset_filtered_combined_archaic_snvs_2.vcf.gz

# index concantenated SNV file
bcftools index -t concat_subset_filtered_combined_archaic_snvs.vcf.gz
echo "indexing complete"
