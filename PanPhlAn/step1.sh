#!/bin/bash
cd /hwfssz5/CNGB_RC/P23Z10200N0829/zhaiqiangrong/fecalmeta/pantrial/Bacteroides_thetaiotaomicron_analysis/

input="fail"




while read -r line;do 
  echo 'export PATH="/hwfssz1/ST_HEALTH/P17Z10200N0306/zhaiqiangrong/minicoda/bin:$PATH" && \' > './shell/bacteroides_thetaiotaomicron_'${line}'_step1.sh'
  echo 'export PATH="/opt/software/conda/envs/mpa/bin/python:$PATH" && \' >> './shell/bacteroides_thetaiotaomicron_'${line}'_step1.sh'
  echo 'export PATH="/opt/software/conda/envs/mpa/bin:$PATH" && \'>> './shell/bacteroides_thetaiotaomicron_'${line}'_step1.sh'
  echo 'export PATH="/opt/software/genomics_tools/samtools/samtools-1.13/bin:$PATH"' >> './shell/bacteroides_thetaiotaomicron_'${line}'_step1.sh'
  
  echo '/opt/software/conda/envs/mpa/bin/python /opt/software/conda/envs/mpa/bin/panphlan_map.py -i /hwfssz5/CNGB_RC/P23Z10200N0829/zhaiqiangrong/fecalmeta/pantrial/Bacteroides_thetaiotaomicron_analysis/FQdir/'${line}'/'${line}'.fastq.bz2 --indexes /hwfssz5/CNGB_RC/P23Z10200N0829/zhaiqiangrong/fecalmeta/pantrial/Bacteroides_thetaiotaomicron/Bacteroides_thetaiotaomicron -p /hwfssz5/CNGB_RC/P23Z10200N0829/zhaiqiangrong/fecalmeta/pantrial/Bacteroides_thetaiotaomicron/Bacteroides_thetaiotaomicron_pangenome.tsv  -o /hwfssz5/CNGB_RC/P23Z10200N0829/zhaiqiangrong/fecalmeta/pantrial/Bacteroides_thetaiotaomicron_analysis/Result/bacteroides_thetaiotaomicron_'${line}'.tsv --tmp /hwfssz5/CNGB_RC/P23Z10200N0829/zhaiqiangrong/fecalmeta/pantrial/TMP/ && \' >> './shell/bacteroides_thetaiotaomicron_'${line}'_step1.sh'
  echo 'rm /hwfssz5/CNGB_RC/P23Z10200N0829/zhaiqiangrong/fecalmeta/pantrial/Bacteroides_thetaiotaomicron_analysis/FQdir/'${line}'/'${line}'.fastq.bz2 ' >> './shell/bacteroides_thetaiotaomicron_'${line}'_step1.sh'

done < "$input"






