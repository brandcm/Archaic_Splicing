#!/bin/bash
#$ -N slop_tabix_grch37_gene_coordinates
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/2_archaic_data_preparation_for_spliceai_analysis/slop_tabix_grch37_gene_coordinates.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/2_archaic_data_preparation_for_spliceai_analysis/slop_tabix_grch37_gene_coordinates.err
#$ -l h_rt=2:00:00
#$ -l mem_free=5G

# load modules
module load CBI
module load bedtools2/2.30.0
module load htslib/1.13

# change directories
cd ../../data/annotations

# perform slop
bedtools slop -i grch37_gene_annotations.bed -g hg19.genome -b 100 > slopped_grch37_gene_annotations.bed
echo "bed file slopped"

# bgzip and tabix bed file
bgzip -c slopped_grch37_gene_annotations.bed > slopped_grch37_gene_annotations.bed.gz
tabix -p bed slopped_grch37_gene_annotations.bed.gz 
echo "bed file bgzipped and indexed"

