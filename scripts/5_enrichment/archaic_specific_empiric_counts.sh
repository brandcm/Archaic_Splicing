#!/bin/bash
#$ -N archaic_specific_empiric_counts
#$ -t 1-5
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/5_enrichment/archaic_specific_empiric_counts.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/5_enrichment/archaic_specific_empiric_counts.err
#$ -l h_rt=24:00:00
#$ -l mem_free=10G

# load modules


# load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate basics

# change directories
cd ../../data/gene_enrichment/observed

# run enrichment
python3 ../../../scripts/5_enrichment/archaic_specific_empiric_counts.py --set_name archaic_specific_2 --iterations 2000 --array_id $SGE_TASK_ID
echo "complete"

# join files
cd ../empiric_counts
join archaic_specific_2_GWAS_empiric_counts_1.txt archaic_specific_2_GWAS_empiric_counts_2.txt > archaic_specific.tmp
for f in archaic_specific_2_GWAS_empiric_counts_3.txt archaic_specific_2_GWAS_empiric_counts_4.txt archaic_specific_2_GWAS_empiric_counts_5.txt; do join archaic_specific.tmp $f > archaic_specific.tmpf && mv archaic_specific.tmpf archaic_specific.tmp; done
mv archaic_specific.tmp archaic_specific_2_GWAS_empiric_counts.txt

join archaic_specific_2_HPO_empiric_counts_1.txt archaic_specific_2_HPO_empiric_counts_2.txt > archaic_specific.tmp
for f in archaic_specific_2_HPO_empiric_counts_3.txt archaic_specific_2_HPO_empiric_counts_4.txt archaic_specific_2_HPO_empiric_counts_5.txt; do join archaic_specific.tmp $f > archaic_specific.tmpf && mv archaic_specific.tmpf archaic_specific.tmp; done
mv archaic_specific.tmp archaic_specific_2_HPO_empiric_counts.txt
