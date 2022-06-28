#!/bin/bash
#$ -N concat_sQTLs
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/analysis/concat_sQTLs.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/analysis/concat_sQTLs.err
#$ -l h_rt=1:00:00
#$ -l mem_free=5G

# load modules


# run script
cd ../../data/GTEx_sQTLs/signif_pairs
for f in *.txt; do awk -F _ '{print $1"_"$2,$3,$4,ARGV[2]}' OFS='\t' $f > chrom_pos_$f & done
wait

for f in chrom_pos_*.txt; do sed 's/\..*//' $f | tail -n +2 > new_$f & done
wait

cat new_*.txt > concat_sQTLs.txt

mv concat_sQTLs.txt ../

rm chrom_pos_*.txt
rm new_*.txt
