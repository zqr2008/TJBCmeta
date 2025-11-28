version 1.0

workflow GenomeQC{
  input{
    File fa_dir
    String file_type="fa"
    String sampleID
    String method
  }
  call run_checkm{
    input:
    fa_dir=fa_dir,
    file_type=file_type,
    sampleID=sampleID,
    method=method
  }
  call run_quast{
    input:
    fa_dir=fa_dir,
    file_type=file_type,
    sampleID=sampleID,
    method=method
  }
  output{
    File checkm_out=run_checkm.checkm_out
    File quast_out=run_quast.quast_out
  }
}
task run_checkm{
  input {
    File fa_dir
    String file_type="fa"
    String sampleID
    String method
  }
  command <<<
    echo $TMPDIR
    export TMPDIR=/tmp
    export PATH=/usr/local/envs/drep_tmp/bin:$PATH

    mkdir ./checkm
    checkm lineage_wf -f ./checkm/output_~{sampleID}_~{method}.txt -t 4 -x ~{file_type} --tab_table ~{fa_dir} ./checkm

  >>>
  runtime {
    docker_url: ""
    req_cpu: 4
    req_memory: "64Gi"
  }
  output {
    File checkm_out = "./checkm"
  }
}

task run_quast{
  input {
    File fa_dir
    String file_type="fa"
    String sampleID
    String method
  }
  command <<<
    echo $TMPDIR
    export TMPDIR=/tmp
    mkdir ./quast
    quast.py ~{fa_dir}/*.~{file_type} -o ./quast/~{sampleID}_~{method}

  >>>
  runtime {
    docker_url: ""
    req_cpu: 4
    req_memory: "16Gi"
  }
  output {
    File quast_out = "./quast"
  }
}