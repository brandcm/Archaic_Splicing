#!/bin/bash
#$ -N spliceai_indel_array
#$ -t 100-201
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/spliceai_indel_array.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/spliceai_indel_array.err
#$ -l h_rt=48:00:00
#$ -l mem_free=20G

# load modules


# load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate spliceai

spliceai -I ../data/split.indel.a$SGE_TASK_ID.vcf.gz -O ../data/split.indel.a$SGE_TASK_ID.vcf -R ../data/hg19.fa -A grch37
