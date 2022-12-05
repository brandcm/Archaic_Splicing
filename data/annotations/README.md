This directory contains annotation files primarily used early in the analysis pipeline; however, a handful are used in later analyses. There are a number of large files missing from this directory such as the GENCODE gff, the gnomAD file, and the phyloP file. I've attempted to note this where relevant but please reach out via email with any questions.

- GENCODE_Release_40_hg38_N_isoforms.txt is a file with the number of predicted isoforms per gene using GENCODE, version 40 for hg38.

- grch37_exon_annotations.txt is the hg19/GRCh37 annotation file from SpliceAI.

- grch37_gene_annotations.bed is a BED file with hg19/GRCh37 gene coordinates using exons defined from GENCODE, version 24.

- grch37_gene_annotations_without_chr.bed is a BED file with hg19/GRCh37 gene coordinates using exons defined from GENCODE, version 24 without the "chr" prefix in the chromosome names.

- grch38_gene_annotations.bed is a BED file with GRCh38 gene coordinates using exons defined from GENCODE, version 24.

- hg19.genome is a chromosome lengths file for hg19/GRCh37 with two tab-delimited fields: 1) the chromosome name and 2) the length in bp.

- hg19_exons.bed.gz is a BED file with hg19 exon coordinates using exons defined from GENCODE, version 24.

- hg19_exons.txt.gz is a file with 1-based hg19 exon coordinates using exons defined from GENCODE, version 24. This file is used to help generate the hg19_start_codons.txt.gz file.

- hg19_start_codons.txt.gz is file with 1-based coordinates for one translation start codon per gene from GENCODE, version 24.

- hg38.genome is a chromosome lengths file for GRCh38 with two tab-delimited fields: 1) the chromosome name and 2) the length in bp.

- major_spliceosome_genes.txt is the output from HGNC for spliceosome associated genes.

- major_spliceosome_genes_hg38.bed is a BED file with hg38 gene coordinates for 140 spliceosome associated genes.

- slopped_grch37_gene_annotations.bed is a BED file with hg19 gene coordinates "buffered" by 100 bp upstream and downstream of the gene. See bedtools "slop".
