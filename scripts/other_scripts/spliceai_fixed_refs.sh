#!/bin/bash
#$ -N spliceai_fixed_refs
#$ -t 1-6
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/other_scripts/spliceai_fixed_refs.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/other_scripts/spliceai_fixed_refs.err
#$ -l h_rt=48:00:00
#$ -l mem_free=30G

# load modules


# load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate spliceai

# change directories
cd ../../data/fixed_refs/vcfs

spliceai -I fixed_reference_"$SGE_TASK_ID".vcf -O fixed_reference_"$SGE_TASK_ID"_annotated.vcf -R ../new_reference/hg19_with_fixed_refs.fa -A grch37
