version 1.0

workflow Binning{
  input{
    String sampleID
    File contig_file
    File contig_fa
    File read1_file
    File read2_file
  }
  
  call Test{
    input:
    sampleID=sampleID
  }
  
  call RunBowtie2{
    input:
    sampleID=sampleID,
    contig_file=contig_file,
    read1_file=read1_file,
    read2_file=read2_file
  }

  call RunMaxbin2{
    input:
    sampleID=sampleID,
    contig_file=contig_file,
    read1_file=read1_file,
    read2_file=read2_file
  }

  call RunMetabat2{
    input:
    sampleID=sampleID,
    contig_fa=contig_fa,
    bowtie2_out=RunBowtie2.bowtie2_out
  }

  call RunConcoct{
    input:
    sampleID=sampleID,
    contig_fa=contig_fa,
    bowtie2_out=RunBowtie2.bowtie2_out
  }

  call RunDASTool{
    input:
    sampleID=sampleID,
    contig_fa=contig_fa,
    maxbin2_out=RunMaxbin2.maxbin2_out,
    metabat2_out=RunMetabat2.metabat2_out,
    concoct_out=RunConcoct.concoct_out
  }

  output{
    File maxbin2_out=RunMaxbin2.maxbin2_out
    File metabat2_out=RunMetabat2.metabat2_out
    File concoct_out=RunConcoct.concoct_out
    File dastool_out=RunDASTool.dastool_out
  }
}

task Test{
  input{
    String sampleID
  }
  command <<<
    echo ~{sampleID} >test.log
  >>>
  runtime {
    docker_url: ""
    req_cpu: 4
    req_memory: "16Gi"
  }
  output {
    File test_out="test.log"
  }
}

task RunBowtie2{
  input {
    String sampleID
    File contig_file
    File read1_file
    File read2_file
  }
  command <<<
    export PATH=/opt/software/conda/envs/dasTool/bin/:/opt/software/conda/envs/maxbin2/bin/:/opt/software/conda/envs/concoct/bin/:/opt/software/conda/envs/metabat2/bin/:$PATH
    mkdir -p ./bowtie2/index ./bowtie2/out
    /opt/software/conda/envs/maxbin2_env/bin/bowtie2-build -f ~{contig_file} ./bowtie2/index/~{sampleID} --threads 4
    /opt/software/conda/envs/maxbin2_env/bin/bowtie2 -1 ~{read1_file} -2 ~{read2_file} -p 4 -x ./bowtie2/index/~{sampleID} -S ./bowtie2/out/~{sampleID}.sam
    /opt/software/conda/envs/concoct/bin/samtools view -@ 4 -b -S ./bowtie2/out/~{sampleID}.sam -o ./bowtie2/out/~{sampleID}.bam
    /opt/software/conda/envs/concoct/bin/samtools sort -@ 4 -l 9 -O BAM ./bowtie2/out/~{sampleID}.bam -o ./bowtie2/out/~{sampleID}.sort.bam
    /opt/software/conda/envs/concoct/bin/samtools index ./bowtie2/out/~{sampleID}.sort.bam
    rm ./bowtie2/out/~{sampleID}.sam ./bowtie2/out/~{sampleID}.bam
  >>>
  runtime {
    docker_url: ""
    req_cpu: 4
    req_memory: "16Gi"
  }
  output {
    File bowtie2_out="./bowtie2"
  }
}

task RunMaxbin2{
  input {
    String sampleID
    File contig_file
    File read1_file
    File read2_file
  }
  command <<<
    export  OMPI_MCA_opal_cuda_support=true
    export  OMPI_MCA_pml="ucx" && export OMPI_MCA_osc="ucx"
    export PATH=/opt/software/conda/envs/maxbin2_env/bin/:/opt/software/programming_language/ActivePerl-5.28.1/bin/:$PATH
    mkdir -p ./maxbin2
    /opt/software/conda/envs/maxbin2_env/bin/run_MaxBin.pl -contig ~{contig_file} -reads ~{read1_file} -reads2 ~{read2_file} -out ./maxbin2/~{sampleID} -preserve_intermediate
    rm ./maxbin2/~{sampleID}.sam0 ./maxbin2/~{sampleID}.sam1
  >>>
  runtime {
    docker_url: ""
    req_cpu: 4
    req_memory: "16Gi"
  }
  output {
    File maxbin2_out="./maxbin2"
  }
}

task RunMetabat2{
  input {
    String sampleID
    File contig_fa
    File bowtie2_out
  }
  command <<<
    export PATH=/opt/software/conda/envs/metabat2/bin/:$PATH
    export LD_LIBRARY_PATH=/opt/software/conda/envs/boost185/lib/:$LD_LIBRARY_PATH
    mkdir -p ./metabat2
    /opt/software/conda/envs/metabat2/bin/jgi_summarize_bam_contig_depths --outputDepth ./metabat2/~{sampleID}.depth.txt ~{bowtie2_out}/out/~{sampleID}.sort.bam
    /opt/software/conda/envs/metabat2/bin/metabat2 -m 1500 -t 4 -i ~{contig_fa} -a ./metabat2/~{sampleID}.depth.txt -o ./metabat2/~{sampleID} -v
  >>>
  runtime {
    docker_url: ""
    req_cpu: 4
    req_memory: "16Gi"
  }
  output {
    File metabat2_out="./metabat2"
  }
}

task RunConcoct{
  input {
    String sampleID
    File contig_fa
    File bowtie2_out
  }
  command <<<
    export PATH=/opt/software/conda/envs/concoct/bin/:$PATH
    mkdir -p ./concoct/out ./concoct/fasta_bins
    /opt/software/conda/envs/concoct/bin/cut_up_fasta.py ~{contig_fa} -c 10000 -o 0 --merge_last -b ./concoct/~{sampleID}_10K.bed > ./concoct/~{sampleID}_10K.fa
    /opt/software/conda/envs/concoct/bin/concoct_coverage_table.py ./concoct/~{sampleID}_10K.bed ~{bowtie2_out}/out/~{sampleID}.sort.bam > ./concoct/~{sampleID}_coverage_table.tsv
    /opt/software/conda/envs/concoct/bin/concoct --composition_file ./concoct/~{sampleID}_10K.fa --coverage_file ./concoct/~{sampleID}_coverage_table.tsv -b ./concoct/out/
    /opt/software/conda/envs/concoct/bin/merge_cutup_clustering.py ./concoct/out/clustering_gt1000.csv > ./concoct/out/clustering_merged.csv
    /opt/software/conda/envs/concoct/bin/extract_fasta_bins.py ~{contig_fa} ./concoct/out/clustering_merged.csv --output_path ./concoct/fasta_bins
  >>>
  runtime {
    docker_url: ""
    req_cpu: 4
    req_memory: "16Gi"
  }
  output {
    File concoct_out="./concoct"
  }
}

task RunDASTool{
  input {
    String sampleID
    File contig_fa
    File maxbin2_out
    File metabat2_out
    File concoct_out
  }
  command <<<
    export PATH=/opt/software/conda/envs/dasTool/bin/:$PATH
    mkdir -p ./DAS_tool
    /opt/software/conda/envs/dasTool/bin/Fasta_to_Contig2Bin.sh -i ~{metabat2_out}/ -e fa | cut -f 1,4 >./DAS_tool/~{sampleID}_metabat2_contigs2bin.tsv
    /opt/software/conda/envs/dasTool/bin/Fasta_to_Contig2Bin.sh -i ~{maxbin2_out}/ -e fasta >./DAS_tool/~{sampleID}_maxbin2_contigs2bin.tsv
    perl -pe "s/,/\tconcoct./g;" ~{concoct_out}/out/clustering_gt1000.csv | sed 's/.concoct_part_[0-9]*//;1d' | sed 's/\t/\t'~{sampleID}_'/' > ./DAS_tool/~{sampleID}_concoct_contigs2bin.tsv
    /opt/software/conda/envs/dasTool/bin/DAS_Tool --search_engine diamond -l concoct,maxbin,metabat -c ~{contig_fa} -o ./DAS_tool/~{sampleID} -i ./DAS_tool/~{sampleID}_concoct_contigs2bin.tsv,./DAS_tool/~{sampleID}_maxbin2_contigs2bin.tsv,./DAS_tool/~{sampleID}_metabat2_contigs2bin.tsv --write_bin_evals --write_bins
  >>>
  runtime {
    docker_url: ""
    req_cpu: 4
    req_memory: "16Gi"
  }
  output {
    File dastool_out="./DAS_tool"
  }
}