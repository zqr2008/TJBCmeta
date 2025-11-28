version 1.0

workflow dRep{
  input{
    Array[File] fa_dir
    String file_type="fa"
  }
  call dRepTask{
    input:
    fa_dir=fa_dir,
    file_type=file_type
  }
  output{
    File result=dRepTask.drep_out
  }
}
task dRepTask{
  input {
    Array[File] fa_dir
    String file_type="fa"
  }
  command <<<
    echo $TMPDIR
    export TMPDIR=/tmp
    export PATH=/usr/local/envs/drep_tmp/bin:$PATH
    mkdir -p ./genome_dir ./drep_out
    /usr/local/envs/drep_tmp/bin/dRep check_dependencies > ./drep_out/check_dependencies.log
    echo "~{sep='\n' fa_dir}" > ./drep_out/dirlist.txt
    cat ./drep_out/dirlist.txt | while read dir; do
      sampID=$(echo ${dir} | sed 's/.*\///;s/_DASTool_bins//' )
      ls ${dir} | while read fa; do 
        cp ${dir}/${fa} ./genome_dir/${sampID}.${fa}
      done
    done
    ls ./genome_dir > ./drep_out/genome_list.log
    /usr/local/envs/drep_tmp/bin/dRep dereplicate -g ./genome_dir/*.~{file_type} -p 8 -l 50000 -comp 50 -con 10 -sa 0.98 --S_algorithm gANI -nc 0.3 -d ./drep_out
  >>>
  runtime {
    docker_url: ""
    req_cpu: 8
    req_memory: "64Gi"
  }
  output {
    File drep_out = "./drep_out"
  }
}
