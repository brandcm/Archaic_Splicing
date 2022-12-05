# change directories
cd ../../data/annotations

# retrieve GFF and format
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/GRCh37_mapping/gencode.v24lift37.basic.annotation.gff3.gz
gzip -d gencode.v24lift37.basic.annotation.gff3.gz
awk -F '[\t,;]' '$3 == "start_codon" {print $1,$2,$4,$5,$7,$15}' OFS='\t' gencode.v24lift37.basic.annotation.gff3 |
awk '{gsub("gene_name=", "", $6); print}' OFS='\t' | awk '!seen[$0]++' | awk '{print $1,$3,$4,$5,$6,$2}' OFS='\t' > tmp.txt

# remove start codons that dont' fall among defined set of exons (in 1-based coordinates)
while IFS= read -r line; do
  echo "$line" > start.txt
  gene=$( echo "$line" | cut -f5 )
  echo "$gene"
  tab=`echo -e "\t"`
  grep "${tab}$gene" hg19_exons.txt > gene_exons.txt
  bedtools intersect -a start.txt -b gene_exons.txt -wa >> intersect.txt; 
done < "tmp.txt"

# let's catch any duplicates
awk '!seen[$0]++' OFS='\t' intersect.txt > tmp2.txt

# forward strand genes
awk '$4 == "+" {print $0}' tmp2.txt | sort -k5,5 -k6,6 -k2n,2 | awk '!seen[$5]++' OFS='\t' > forward_tmp.txt

# reverse strand genes
awk '$4 == "-" {print $0}' tmp2.txt | sort -k5,5 -k6,6 -k3nr,3 | awk '!seen[$5]++' OFS='\t' > reverse_tmp.txt

# cat and sort
cat forward_tmp.txt reverse_tmp.txt > tmp3.txt
sort -k1,1 -k2n,2 tmp3.txt | awk '{print $1,$6,$2,$3,$4,$5}' OFS='\t' > hg19_start_codons.txt

# remove temp files
rm tmp.txt && rm start.txt && rm gene_exons.txt && rm intersect.txt && rm tmp2.txt && rm forward_tmp.txt && rm reverse_tmp.txt && rm tmp3.txt
