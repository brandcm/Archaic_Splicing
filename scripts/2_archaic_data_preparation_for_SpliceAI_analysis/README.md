This directory houses scripts to futher prepare the archaic data for SpliceAI analysis. Scripts are listed in the order in which they should be launched.

- slop_tabix_grch37_gene_coordinates.sh generates a BED file from which to subset archaic variants. SpliceAI does not annotate variants that lie outside of gene bodies. Thus, the resulting BED file reflects the gene start and stop coordinates with a 100 bp "buffer" on either end.

- subset_vcfs.sh combines the per contig SNV files and subsets the filtered data to those variants located within or nearby (+- 100 bp) gene bodies, as annoated by GENCODE, Human Release 24 for hg19.

- split_*_vcf.sh splits the InDels and SNVs into smaller files for SpliceAI to more efficiently process. The current split is set at 5,000 variants but this can be adjusted. 
