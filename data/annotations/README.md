This directory contains annotation files primarily used early in the analysis pipeline; however, a handful are used in later analyses. There are a number of large files missing from this directory such as the GENCODE gff, the gnomAD file, and the phyloP file. I've attempted to note this where relevant but please reach out via email with any questions.

- hg19.genome is a chromosome lengths file for hg19/GRCh37 with two tab-delimited fields: 1) the chromosome name and 2) the length in bp.

- hg19_exons.bed.gz is a BED file with hg19 exon coordinates using exons defined from GENCODE, version 24.

- hg19_start_codons.txt.gz is file with 1-based coordinates for one translation start codon per gene from GENCODE, version 24.

- hg38.genome is a chromosome lengths file for GRCh38 with two tab-delimited fields: 1) the chromosome name and 2) the length in bp.

- major_spliceosome_genes.txt is the output from HGNC for spliceosome associated genes.

- major_spliceosome_genes_hg38.bed is a BED file with hg38 gene coordinates for 140 spliceosome associated genes.

- slopped_grch37_gene_annotations.bed is a BED file with hg19 gene coordinates "buffered" by 100 bp upstream and downstream of the gene. See bedtools "slop".
