#!/bin/bash
#$ -N spliceai_introgressed_refs
#$ -t 1-2
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/other_scripts/spliceai_introgressed_refs.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/other_scripts/spliceai_introgressed_refs.err
#$ -l h_rt=48:00:00
#$ -l mem_free=20G

# load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate spliceai

# change directories
cd ../../data/introgressed_refs/vcfs

spliceai -I introgressed_refs_"$SGE_TASK_ID".vcf -O introgressed_refs_"$SGE_TASK_ID"_annotated.vcf -R ../new_reference/hg19_with_introgressed_alts.fa -A grch37
