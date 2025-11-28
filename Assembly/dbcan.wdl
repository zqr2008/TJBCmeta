version 1.0

workflow CAZY{
  input{
    File drep_prodigal_dir
    String subjectID
  }
  
  call CAZy{
    input:
    drep_prodigal_dir=drep_prodigal_dir,
    subjectID=subjectID
  }

  output{
    File out=CAZy.dbcan_out
  }
}

task CAZy{
  input {
    File drep_prodigal_dir
    String subjectID
  }
  command <<<

    source /opt/software/conda/envs/dbcan/bin/activate
    ln -s /opt/DB/dbCAZy ./db

    ls ~{drep_prodigal_dir}/*.faa | while read faa; do
      genome=$(basename ${faa} | sed 's/.fa.faa//')
      mkdir ./${genome}.dbcan_out
      run_dbcan ${faa} protein --out_dir ./${genome}.dbcan_out
      #run_dbcan ${fa} prok --out_dir ./${genome}.dbcan_out
    done

    tar zcvf ./~{subjectID}.dbcan_out.tar.gz ./*.dbcan_out

    source /opt/software/conda/envs/dbcan/bin/deactivate

  >>>
  runtime {
    docker_url: ""
    req_cpu: 1
    req_memory: "8Gi"
  }
  output {
    File dbcan_out = "./~{subjectID}.dbcan_out.tar.gz"
  }
}