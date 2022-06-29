This directory contains scripts to apply SpliceAI to archaic variants and concatenate the results.

- concat_spliceai_*_vcfs.sh concatenates the split outputs for each set. InDels and SNVs are left separate and are later joined in the [Generate Archaic Dataframe notebook](https://github.com/brandcm/Archaic_Splicing/blob/main/scripts/notebooks/1_generate_archaic_dataframe.ipynb). Also note that I split variants with multiple annotations as the deltas do not necessarily agree across annotations (i.e., one variant may have no splicing effect on gene A yet affect gene B). This results in a dataframe with the same number of fields. Therefore, the length of the resulting dataframe(s) does not equal the number of positions represented in the file(s).

- spliceai_*_array.sh launches an array job in which SpliceAI is applied to either the archaic InDels or SNVs. Note that the array size depends on the number of split files.
