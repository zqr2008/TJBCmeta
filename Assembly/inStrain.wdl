version 1.0

workflow AssemblyProfiling{
  input{
    File drep_genomes_dir
    File drep_prodigal_dir
    String genome_type="fa"
    Array[String] samples
    Array[File] read1_list
    Array[File] read2_list
    String subjectID
  }

  call checkSampleInfo{
    input:
    samples=samples,
    read1_list=read1_list,
    read2_list=read2_list
  }
  
  call processFA{
    input:
    drep_genomes_dir=drep_genomes_dir,
    drep_prodigal_dir=drep_prodigal_dir,
    genome_type=genome_type
  }

  call prodigalCat{
    input:
    drep_genomes_dir=processFA.genomes_dir,
    drep_prodigal_dir=processFA.prodigal_dir,
    subjectID=subjectID,
    genome_type=genome_type
  }

  call fastaCat{
    input:
    drep_genomes_dir=processFA.genomes_dir,
    subjectID=subjectID,
    genome_type=genome_type
  }
  
  call parse_stb{
    input:
    drep_genomes_dir=processFA.genomes_dir,
    subjectID=subjectID,
    genome_type=genome_type
  }
  
  call build_bowtie2_index{
    input:
    fasta_file=fastaCat.MultiFasta,
    subjectID=subjectID
  }

  scatter (i in range(length(samples))) {
    call RunBowtie2{
      input:
      ref_index=build_bowtie2_index.ref_index,
      subjectID=subjectID,
      read1_file=read1_list[i],
      read2_file=read2_list[i],
      sampleID=samples[i]
    }
    call inStrainProfile{
      input:
      bam_file=RunBowtie2.bam_file,
      bai_file=RunBowtie2.bai_file,
      read1_file=read1_list[i],
      read2_file=read2_list[i],
      sampleID=samples[i],
      subjectID=subjectID,
      MultiFasta=fastaCat.MultiFasta,
      prodigalMultiFasta=prodigalCat.prodigalMultiFasta,
      stb_file=parse_stb.stb_file
    }
  }

  if (length(samples) > 1){
    call inStrainCompare{
      input:
      bam_list=RunBowtie2.bam_file,
      bai_list=RunBowtie2.bai_file,
      profile_list=inStrainProfile.inStrain_out,
      subjectID=subjectID,
      stb_file=parse_stb.stb_file
    }
  }
  
  output{
    File genome_list=processFA.genome_list
    File prodigal_list=processFA.prodigal_list
    Array[File] inStrain_profile=inStrainProfile.inStrain_out_out
    File prodigalMultiFasta=prodigalCat.prodigalMultiFasta
    File MultiFasta=fastaCat.MultiFasta
    File stb_file=parse_stb.stb_file
    File? inStrain_compare=inStrainCompare.inStrain_compare
    File check_sample_info=checkSampleInfo.check_sample_info
    Array[File] checkFile=RunBowtie2.checkFile
  }
}

task checkSampleInfo{
  input {
    Array[String] samples
    Array[File] read1_list
    Array[File] read2_list
  }
  command <<<
    echo "~{sep='\n' samples}" >samples.txt
    echo "~{sep='\n' read1_list}" >read1_list.txt
    echo "~{sep='\n' read2_list}" >read2_list.txt

    paste samples.txt read1_list.txt read2_list.txt > ./check_sample_info.txt
  >>>
  runtime {
    docker_url: ""
    req_cpu: 1
    req_memory: "1Gi"
  }
  output {
    File check_sample_info = "./check_sample_info.txt"
  }
}

task processFA{
  input {
    File drep_genomes_dir
    File drep_prodigal_dir
    String genome_type="fa"
  }
  command <<<
    mkdir -p ./dereplicated_genomes
    ls ~{drep_genomes_dir}/*.~{genome_type} | while read fa; do 
      sampleID=$(basename ${fa} | sed 's/dereplicated_genomes\.//' | awk -F "." '{print $1}')
      sed 's/k77/'${sampleID}'_k77/' ${fa} >./dereplicated_genomes/$(basename ${fa})
    done
    ls ./dereplicated_genomes > dereplicated_genomes.list
    
    mkdir -p ./prodigal
    ls dereplicated_genomes/*.~{genome_type} | while read fa; do 
      sampleID=$(basename ${fa} | sed 's/dereplicated_genomes\.//' | awk -F "." '{print $1}')
      sed 's/k77/'${sampleID}'_k77/' ~{drep_prodigal_dir}/$(basename ${fa}).fna >./prodigal/$(basename ${fa}).fna
    done
    ls ./prodigal > prodigal.list
  >>>
  runtime {
    docker_url: ""
    req_cpu: 1
    req_memory: "1Gi"
  }
  output {
    File genome_list="dereplicated_genomes.list"
    File prodigal_list="prodigal.list"
    File genomes_dir = "./dereplicated_genomes"
    File prodigal_dir = "./prodigal"
  }
}

task prodigalCat{
  input {
    File drep_genomes_dir
    File drep_prodigal_dir
    String subjectID
    String genome_type="fa"
  }
  command <<<
    rm -f ./~{subjectID}.prodigal.multi.fna
    ls ~{drep_genomes_dir}/*.~{genome_type} > ./drep_genomes_list.txt
    cat drep_genomes_list.txt | sed 's/.*\///' | while read fa; do
      cat ~{drep_prodigal_dir}/${fa}.fna >> ./~{subjectID}.prodigal.multi.fna
    done
  >>>
  runtime {
    docker_url: ""
    req_cpu: 1
    req_memory: "1Gi"
  }
  output {
    File prodigalMultiFasta = "./~{subjectID}.prodigal.multi.fna"
    File drep_genomes_list = "./drep_genomes_list.txt"
  }
}

task fastaCat{
  input {
    File drep_genomes_dir
    String subjectID
    String genome_type="fa"
  }
  command <<<
    rm -f ./~{subjectID}.multi.fa
    ls ~{drep_genomes_dir}/*.~{genome_type} > ./drep_genomes_list.txt
    cat drep_genomes_list.txt | while read fa; do
      cat ${fa} >> ./~{subjectID}.multi.fa
    done
  >>>
  runtime {
    docker_url: ""
    req_cpu: 1
    req_memory: "1Gi"
  }
  output {
    File MultiFasta = "./~{subjectID}.multi.fa"
    File drep_genomes_list = "./drep_genomes_list.txt"
  }
}

task parse_stb{
  input {
    File drep_genomes_dir
    String subjectID
    String genome_type="fa"
  }
  command <<<
    echo $TMPDIR
    export TMPDIR=/tmp
    export PATH=/usr/local/envs/drep_tmp/bin:$PATH
    parse_stb.py --reverse -f ~{drep_genomes_dir}/*.~{genome_type} -o ./~{subjectID}_cat_genomes.stb
  >>>
  runtime {
    docker_url: ""
    req_cpu: 1
    req_memory: "16Gi"
  }
  output {
    File stb_file = "./~{subjectID}_cat_genomes.stb"
  }
}


task build_bowtie2_index{
  input {
    File fasta_file
    String subjectID
  }
  command <<<
    export PATH=/opt/software/conda/envs/dasTool/bin/:/opt/software/conda/envs/maxbin2/bin/:/opt/software/conda/envs/concoct/bin/:/opt/software/conda/envs/metabat2/bin/:$PATH
    mkdir -p ./ref_index
    /opt/software/conda/envs/maxbin2_env/bin/bowtie2-build --threads 4 ~{fasta_file} ./ref_index/~{subjectID}_ref
  >>>
  runtime {
    docker_url: ""
    req_cpu: 4
    req_memory: "16Gi"
  }
  output {
    File ref_index = "./ref_index"
  }
}

task RunBowtie2{
  input {
    File ref_index
    String subjectID
    File read1_file
    File read2_file
    String sampleID
  }
  command <<<
    echo -e "~{sampleID}\t~{read1_file}\t~{read2_file}" >~{sampleID}.bowtie2_check.txt
    
    export PATH=/opt/software/conda/envs/dasTool/bin/:/opt/software/conda/envs/maxbin2/bin/:/opt/software/conda/envs/concoct/bin/:/opt/software/conda/envs/metabat2/bin/:$PATH
    /opt/software/conda/envs/maxbin2_env/bin/bowtie2 -1 ~{read1_file} -2 ~{read2_file} -p 4 -x ~{ref_index}/~{subjectID}_ref -S ./~{subjectID}.~{sampleID}.sam
    /opt/software/conda/envs/concoct/bin/samtools view -@ 4 -b -S ./~{subjectID}.~{sampleID}.sam -o ./~{subjectID}.~{sampleID}.bam
    /opt/software/conda/envs/concoct/bin/samtools sort -@ 4 -l 9 -O BAM ./~{subjectID}.~{sampleID}.bam -o ./~{subjectID}.~{sampleID}.sort.bam
    /opt/software/conda/envs/concoct/bin/samtools index ./~{subjectID}.~{sampleID}.sort.bam
    rm ./~{subjectID}.~{sampleID}.sam ./~{subjectID}.~{sampleID}.bam
  >>>
  runtime {
    docker_url: ""
    req_cpu: 4
    req_memory: "16Gi"
  }
  output {
    File bam_file = "./~{subjectID}.~{sampleID}.sort.bam"
    File bai_file = "./~{subjectID}.~{sampleID}.sort.bam.bai"
    File checkFile = "~{sampleID}.bowtie2_check.txt"
  }
}

task inStrainProfile{
  input {
    File bam_file
    File bai_file
    File read1_file
    File read2_file
    String sampleID
    String subjectID
    File MultiFasta
    File prodigalMultiFasta
    File stb_file
  }
  command <<<
    echo $TMPDIR
    export TMPDIR=/tmp
    export PATH=/usr/local/envs/inStrain/bin/:$PATH
    mkdir ./bamdir
    ln -s ~{bam_file} ./bamdir/~{subjectID}.~{sampleID}.sort.bam
    ln -s ~{bai_file} ./bamdir/~{subjectID}.~{sampleID}.sort.bam.bai

    inStrain profile ./bamdir/~{subjectID}.~{sampleID}.sort.bam ~{MultiFasta} -o ./~{subjectID}_~{sampleID}.IS -p 4 -g ~{prodigalMultiFasta} -s ~{stb_file}
    
    cp -r ./~{subjectID}_~{sampleID}.IS ./~{subjectID}_~{sampleID}.IS_out
    rm -rf ./~{subjectID}_~{sampleID}.IS_out/raw_data
    
    if [ -e ./~{subjectID}_~{sampleID}.IS/output/~{subjectID}_~{sampleID}.IS_genome_info.tsv ]; then echo done; else exit 1; fi
  >>>
  runtime {
    docker_url: ""
    req_cpu: 4
    req_memory: "112Gi"
  }
  output {
    File inStrain_out="./~{subjectID}_${sampleID}.IS"
    File inStrain_out_out="./~{subjectID}_~{sampleID}.IS_out"
  }
}

task inStrainCompare{
  input {
    Array[File] bam_list
    Array[File] bai_list
    Array[File] profile_list
    String subjectID
    File stb_file
  }
  command <<<
    echo $TMPDIR
    export TMPDIR=/tmp
    export PATH=/usr/local/envs/inStrain/bin/:$PATH
    mkdir ./bamdir ./profiles

    echo "~{sep='\n' bam_list}" >bam_list.txt
    echo "~{sep='\n' bai_list}" >bai_list.txt
    echo "~{sep='\n' profile_list}" >profile_list.txt
    
    cat bam_list.txt | while read bam; do ln -s ${bam} ./bamdir/; done
    cat bai_list.txt | while read bai; do ln -s ${bai} ./bamdir/; done
    cat profile_list.txt | while read profile; do ln -s ${profile} ./profiles/; done

    inStrain compare -i ./profiles/*.IS/ -bams ./bamdir/*.sort.bam -o ./~{subjectID}.IS.COMPARE -p 8 --store_mismatch_locations -s ~{stb_file} --database_mode
  >>>
  runtime {
    docker_url: ""
    req_cpu: 8
    req_memory: "112Gi"
  }
  output {
    File inStrain_compare="./~{subjectID}.IS.COMPARE"
  }
}