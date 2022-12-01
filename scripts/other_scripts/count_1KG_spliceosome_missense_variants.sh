# retrieve information from 1KG .vep output files from Ensembl's Variant Effect Predictor

inds=('HG00137' 'HG00338' 'HG00619' 'HG01198' 'HG01281' 'HG01524' 'HG01802' 'HG02142' 'HG02345' 'HG02629' 'HG03060' 'HG03190' 'HG03708' 'HG03711' 'HG03800' 'HG04014' 'NA11830' 'NA18552' 'NA18868' 'NA19011' 'NA19452' 'NA19741' 'NA20537' 'NA21141' )

cd ../../data/thousand_genomes/spliceosome_variants

echo "starting missense subset"
for ind in ${inds[@]}; do echo "$ind" & awk '($7=="missense_variant") {print $0}' "$ind"_spliceosome_variants.vep > "$ind"_missense_spliceosome_variants.vep; done

echo "starting transcript count"
for ind in ${inds[@]}; do echo "$ind" & wc -l "$ind"_missense_spliceosome_variants.vep; done 

echo "starting variant count"  
for ind in ${inds[@]}; do echo "$ind" & awk '!_[$2]++' "$ind"_missense_spliceosome_variants.vep | wc -l; done

echo "starting variants with a damaging PolyPhen score count" 
for ind in ${inds[@]}; do echo "$ind" & grep 'damaging' "$ind"_missense_spliceosome_variants.vep | awk '!_[$2]++' | wc -l; done

echo "starting variants with a deleterious SIFT score count"
for ind in ${inds[@]}; do echo "$ind" & grep 'deleterious' "$ind"_missense_spliceosome_variants.vep | awk '!_[$2]++' | wc -l; done

echo "starting variants with a damaging PolyPhen score AND a deleterious SIFT score count"
for ind in ${inds[@]}; do echo "$ind" & grep 'damaging' "$ind"_missense_spliceosome_variants.vep | grep 'deleterious' | awk '!_[$2]++' | wc -l; done