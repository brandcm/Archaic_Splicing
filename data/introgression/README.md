This directory contains two introgression callsets from Browning et al. 2018 and Vernot et al. 2016.

- Browning_et_al_2018_Neanderthal_introgressed_variants.txt.gz is a concatenation of all non-African population introgression calls from Browning et al. 2018 (retrieved from https://data.mendeley.com/datasets/y7hyt83vxr/1). I retained only rows with a matching Neanderthal allele.

- Vernot_et_al_2016_introgressed_tag_snps.bed.gz is a concatenation of the ASN, EUR, PNG, and SAS tag SNP files from Vernot et al. 2016 (retrieved from https://drive.google.com/drive/folders/0B9Pc7_zItMCVM05rUmhDc0hkWmc?resourcekey=0-zwKyJGRuooD9bWPRZ0vBzQ) generated using [this script](https://github.com/brandcm/Archaic_Splicing/blob/main/scripts/4_modern_data_preparation/concat_introgressed_variants.sh). Duplicate entries that resulted from merging data for multiple populations were removed. We used the all_tag_snps.*.merged.ALL.0.3_R2_cluster.1KG_phase3_essentials.bed set of files.
