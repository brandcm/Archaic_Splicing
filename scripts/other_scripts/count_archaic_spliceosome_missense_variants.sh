# retrieve information from archaic .vep output files from Ensembl's Variant Effect Predictor

inds=('altai' 'chagyrskaya' 'denisovan' 'vindija' )

cd ../../data/spliceosome

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