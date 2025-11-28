version 1.0

workflow Megahit{
  input{
    File rmhost1
    File rmhost2
    String sampleid
  }
  call megahit{
    input:
    rmhost1=rmhost1,
    rmhost2=rmhost2,
    sampleid=sampleid
  }
  output{
    File fastg = megahit.fastg
    File fagz = megahit.fagz
    File fa = megahit.fa
  }
}


task megahit {
    input {
        File rmhost1
        File rmhost2
        String sampleid
    }
    command {
        
        /opt/software/conda/envs/megahit/bin/megahit -1 ${rmhost1} -2 ${rmhost2} -t 8 --k-list 21,33,55,77 --min-contig-len 100 --out-dir ./"${sampleid}.megahit.out" --out-prefix ${sampleid} 2> "${sampleid}.megahit.log"
        /opt/software/conda/envs/megahit/bin/megahit_toolkit contig2fastg 77 ./"${sampleid}.megahit.out"/"${sampleid}.contigs.fa" > ./"${sampleid}.megahit.out"/"${sampleid}.megahit.scaftigs.fastg"
        /opt/software/conda/envs/megahit/bin/pigz -p 8 "./${sampleid}.megahit.out/${sampleid}.megahit.scaftigs.fastg"
        /opt/software/conda/envs/megahit/bin/pigz -p 8 -k "./${sampleid}.megahit.out/${sampleid}.contigs.fa"
    }
    output {
        File fastg = "./${sampleid}.megahit.out/${sampleid}.megahit.scaftigs.fastg.gz"
        File fagz = "./${sampleid}.megahit.out/${sampleid}.contigs.fa.gz"
        File fa = "./${sampleid}.megahit.out/${sampleid}.contigs.fa"
  }
    runtime {
        docker_url: ""
        req_cpu:8
        req_memory:"64Gi"
  }
}