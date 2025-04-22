#!/bin/bash
set -euo pipefail
echo ==========start at : `date` ==========

# print the potential clades 
#strainphlan -d ./database/mpa_vOct22_CHOCOPhlAnSGB_202212.pkl -s *pkl -o ./db_marker/ --print_clades_only > ./samples.marker.log
#grep "t__" ./samples.marker.log |cut -f4 -d':'|sed 's/\t//g' > ./samples.marker.list

# generate scripts for transmission determination
cat samples.marker.list | while read t; do 

cat > ./shell/${t}.accurate.sh <<EOF
#!bin/bash
mkdir -p ./db_marker/${t} ./Accurateoutput/${t}

#Assessment of person-to-person microbiome strain sharing
#Reference: https://github.com/biobakery/MetaPhlAn/wiki/Strain-Sharing-Inference

#Extracting the species marker genes
extract_markers.py -d ./database/mpa_vOct22_CHOCOPhlAnSGB_202212.pkl -c ${t} -o ./db_marker/${t}

#Strain-level profiling
strainphlan -d ./database/mpa_vOct22_CHOCOPhlAnSGB_202212.pkl -s ./Sample/*pkl -m ./db_marker/${t}/${t}.fna -o ./Accurateoutput/${t} -n 2 -c ${t} --mutation_rates --marker_in_n_samples 10 --sample_with_n_markers 10 --secondary_sample_with_n_markers 10 --phylophlan_mode accurate

#Extraction of pairwise phylogenetic distances
tree_pairwisedists.py -n ./Accurateoutput/${t}/RAxML_bestTree.${t}.StrainPhlAn4.tre ./Accurateoutput/${t}/${t}_nGD.tsv

#Getting species-specific strain identity thresholds
Rscript ./threshold_calc.R ./1w8.metadata.forR ./Accurateoutput/${t}/${t}_nGD.tsv ./Accurateoutput/${t}/

#Detection of strain sharing events
thre=\$(cat ./Accurateoutput/${t}/threshold.txt)
strain_transmission.py --tree ./Accurateoutput/${t}/RAxML_bestTree.${t}.StrainPhlAn4.tre --metadata ./1w8.metadata.forTrans --output_dir ./Accurateoutput/${t} --threshold \${thre}

EOF

cat > ./shell/${t}.fast.sh <<EOF
#!bin/bash
mkdir -p ./Fastoutput/

#use species marker genes extracted before
#Strain-level profiling with "fast" mode
strainphlan -d ./database/mpa_vOct22_CHOCOPhlAnSGB_202212.pkl -s ./Sample/*pkl -m ./db_marker/${t}/${t}.fna -o ./Fastoutput/${t} -n 2 -c ${t} --mutation_rates --marker_in_n_samples 10 --sample_with_n_markers 10 --secondary_sample_with_n_markers 10 --phylophlan_mode fast

#Extraction of pairwise phylogenetic distances
tree_pairwisedists.py -n ./Fastoutput/${t}/RAxML_bestTree.${t}.StrainPhlAn4.tre ./Fastoutput/${t}/${t}_nGD.tsv

#Getting species-specific strain identity thresholds
Rscript ./threshold_calc.R ./1w8.metadata.forR ./Fastoutput/${t}/${t}_nGD.tsv ./Fastoutput/${t}/

#Detection of strain sharing events
thre=\$(cat ./Fastoutput/${t}/threshold.txt)
strain_transmission.py --tree ./Fastoutput/${t}/RAxML_bestTree.${t}.StrainPhlAn4.tre --metadata ./1w8.metadata.forTrans --output_dir ./Fastoutput/${t} --threshold \${thre}

EOF

done

echo ==========end at : `date` ==========
