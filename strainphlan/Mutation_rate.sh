#!bin/bash

s=s__Phocaeicola_vulgatus
t=t__SGB1814
awk -v OFS="\t" 'NR==1{for(i=2;i<=NF;i++){a[i]=$i}}NR!=1{for(j=2;j<=NF;j++){print $1,a[j],$j}}' ./output/${t}.mutation | grep -v "/" >./output/${t}.mutation.reformat
awk -v OFS="\t" 'NR==FNR{s[$1]=$3"\t"$5; f[$1]=$4}NR!=FNR{print $1,s[$1],f[$1],$2,s[$2],f[$2],$3}' ../1w8.metadata.forR ./output/${t}.mutation.reformat | awk '{if($2==$6 && $3=="mother"){print "'$s'\t'${t}'\t"$0"\tIntra-individual(mother)"} else if($2==$6 && $3=="offspring"){print "'$s'\t'${t}'\t"$0"\tIntra-individual(offspring)"}else if($2!=$6 && $4==$8 && ($3"-"$7=="mother-offspring" || $3"-"$7=="offspring-mother")){print "'$s'\t'${t}'\t"$0"\tIntra-family"}else if ($4!=$8 && ($3"-"$7=="mother-offspring" || $3"-"$7=="offspring-mother")){print "'$s'\t'${t}'\t"$0"\tInter-family"}else{if($3"-"$7 == "mother-mother"){print "'$s'\t'${t}'\t"$0"\tInter-mother"}else if($3"-"$7 == "offspring-offspring"){"'$s'\t'${t}'\t"$0"\tInter-offspring"}else{print "'$s'\t'${t}'\t"$0"\tNone"}}}' > ./output/${t}.mutation.reformat.addMeta
awk -v OFS="\t" 'NR==FNR{info[$1]=$2"\t"$3}NR!=FNR{if($3 in info && $7 in info){print $1,$2,$3,$4,$5,$6,info[$3],$7,$8,$9,$10,info[$7],$11,$12}}' pheno.tsv ./output/${t}.mutation.reformat.addMeta > ./output/${t}.mutation.reformat.addMeta.addPheno

grep "Intra-family" ./output/${t}.mutation.reformat.addMeta.addPheno > ./output/${t}.mutation.reformat.addMeta.addPheno.intraFamily

awk 'BEGIN{print "Species\tStrain\tSample1\tSubject1\tRole1\tFamily1\tVisit1\tTraj1\tSample2\tSubject2\tRole2\tFamily2\tVisit2\tTraj2\tMutation_Rate\tGroup\tVisit_pair"}{if($5"-"$11=="mother-offspring"){print $0"\t"$7"-"$13}else if($5"-"$11=="offspring-mother"){print $0"\t"$13"-"$7}else{print $0"\tError"}}' ./output/${t}.mutation.reformat.addMeta.addPheno.intraFamily > ./output/${t}.mutation.reformat.addMeta.addPheno.intraFamily.addVisGroup