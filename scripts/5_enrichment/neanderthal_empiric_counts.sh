#!/bin/bash
#$ -N neanderthal_empiric_counts
#$ -t 2
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/5_enrichment/neanderthal_empiric_counts.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/5_enrichment/neanderthal_empiric_counts.err
#$ -l h_rt=24:00:00
#$ -l mem_free=10G

# load modules


# load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate basics

# change directories
cd ../../data/gene_enrichment/observed

# run enrichment
python3 ../../../scripts/5_enrichment/neanderthal_empiric_counts.py --set_name neanderthal_2 --iterations 2000 --array_id $SGE_TASK_ID
echo "complete"

# join files
cd ../empiric_counts
join neanderthal_2_GWAS_empiric_counts_1.txt neanderthal_2_GWAS_empiric_counts_2.txt > neanderthal.tmp
for f in neanderthal_2_GWAS_empiric_counts_3.txt neanderthal_2_GWAS_empiric_counts_4.txt neanderthal_2_GWAS_empiric_counts_5.txt; do join neanderthal.tmp $f > neanderthal.tmpf && mv neanderthal.tmpf neanderthal.tmp; done
mv neanderthal.tmp neanderthal_2_GWAS_empiric_counts.txt

join neanderthal_2_HPO_empiric_counts_1.txt neanderthal_2_HPO_empiric_counts_2.txt > neanderthal.tmp
for f in neanderthal_2_HPO_empiric_counts_3.txt neanderthal_2_HPO_empiric_counts_4.txt neanderthal_2_HPO_empiric_counts_5.txt; do join neanderthal.tmp $f > neanderthal.tmpf && mv neanderthal.tmpf neanderthal.tmp; done
mv neanderthal.tmp neanderthal_2_HPO_empiric_counts.txt