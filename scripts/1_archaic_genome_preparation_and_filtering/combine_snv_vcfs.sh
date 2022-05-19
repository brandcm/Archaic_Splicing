#!/bin/bash
#$ -N combine_snv_vcfs
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/1_archaic_genome_preparation_and_filtering/combine_snv_vcfs.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/1_archaic_genome_preparation_and_filtering/combine_snv_vcfs.err
#$ -l h_rt=10:00:00
#$ -l scratch=50G

# load modules
module load CBI
module load bcftools/1.13

# copy needed files into the temporary directory
cd "$TMPDIR"
cp ~/../../../group/capra/data/altai_neanderthal/2021-11-17/*.vcf.gz . &
cp ~/../../../group/capra/data/altai_neanderthal/2021-11-17/*.vcf.gz.tbi . &
cp ~/../../../group/capra/data/chagyrskaya_neanderthal/2021-11-17/*.vcf.gz . &
cp ~/../../../group/capra/data/chagyrskaya_neanderthal/2021-11-17/*.vcf.gz.tbi . &
cp ~/../../../group/capra/data/denisovan/2021-11-18/*.vcf.gz . &
cp ~/../../../group/capra/data/denisovan/2021-11-18/*.vcf.gz.tbi . &
cp ~/../../../group/capra/data/vindija_neanderthal/2021-11-18/*.vcf.gz . &
cp ~/../../../group/capra/data/vindija_neanderthal/2021-11-18/*.vcf.gz.tbi . &
cp ~/../../../group/capra/projects/archaic_splicing/data/archaic_genomes/new_chr_names.txt . &
cp ~/../../../group/capra/projects/archaic_splicing/data/archaic_genomes/new_sample_names.txt .
wait
touch *.tbi
echo "moving files finished"

# list of chromosomes
chrs=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22' 'chrX' )

# chagyrskaya does not appear to be properly bgzipped; let's change that
for c in ${chrs[@]}; do bcftools view chagyrskaya_$c.vcf.gz -O z -o bgzipped_chagyrskaya_$c.vcf.gz & done
wait
echo "chagyrskaya bgzipped"

for c in ${chrs[@]}; do bcftools index -t bgzipped_chagyrskaya_$c.vcf.gz & done
wait
echo "indexing chagryskaya done"

# combine VCFs
for c in ${chrs[@]}; do bcftools merge altai_$c.vcf.gz bgzipped_chagyrskaya_$c.vcf.gz denisovan_$c.vcf.gz vindija_$c.vcf.gz -O z -o combined_$c.vcf.gz & done
wait
echo "combining VCFs complete"

# rename chromosomes and samples, index new VCFs
for c in ${chrs[@]}; do bcftools annotate --rename-chrs new_chr_names.txt combined_$c.vcf.gz -o combined_2_$c.vcf.gz & done
wait
echo "chromosomes renamed"

for c in ${chrs[@]}; do bcftools reheader -s new_sample_names.txt combined_2_$c.vcf.gz -o combined_archaic_$c.vcf.gz & done
wait
echo "samples renamed"

for c in ${chrs[@]}; do bcftools index -t combined_archaic_$c.vcf.gz & done
wait
echo "combined VCFs indexed"

# move new VCFs and indices
mv combined_archaic_*.vcf.gz ~/../../../group/capra/projects/archaic_splicing/data/archaic_genomes
mv combined_archaic_*.vcf.gz.tbi ~/../../../group/capra/projects/archaic_splicing/data/archaic_genomes
echo "script completed"
