#!/bin/bash
#$ -N spliceai_snv_array
#$ -t 100-443
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/3_spliceai_analysis/spliceai_snv_array.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/3_spliceai_analysis/spliceai_snv_array.err
#$ -l h_rt=48:00:00
#$ -l mem_free=20G

# load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate spliceai

# set path to hg19 and change directories
hg19="../../../../data/hg19_fasta/2022-03-14/hg19.fa"
cd ../../data/archaic_genomes

spliceai -I split.x"$SGE_TASK_ID".vcf.gz -O ../spliceai_outputs/split.x"$SGE_TASK_ID".vcf -R "$hg19" -A grch37
