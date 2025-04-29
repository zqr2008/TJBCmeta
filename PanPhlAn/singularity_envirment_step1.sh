#!/bin/bash
cd /hwfssz5/CNGB_RC/P23Z10200N0829/zhaiqiangrong/fecalmeta/pantrial/Bacteroides_thetaiotaomicron_analysis/

input="fail"




while read -r line;do 
  
  echo '/share/app/singularity/3.8.1/bin/singularity exec --cleanenv --bind /hwfssz1,/zfssz8,/jdfssz1,/hwfssz5/ /hwfssz5/CNGB_RC/P23Z10200N0829/zhaiqiangrong/fecalmeta/pantrial/siffile/dockerhub.centos7.9.2009.development.basic.omicsbasic.meta_v2.1.sif bash /hwfssz5/CNGB_RC/P23Z10200N0829/zhaiqiangrong/fecalmeta/pantrial/Bacteroides_thetaiotaomicron_analysis/shell/bacteroides_thetaiotaomicron_'${line}'_step1.sh ' > './shell/bacteroides_thetaiotaomicron_'${line}'_sing.sh'


done < "$input"






