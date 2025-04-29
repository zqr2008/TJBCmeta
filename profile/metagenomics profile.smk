##############################################################
### chenjunhong@genomics.cn
### 2023/04/06
### metagenomics profile v1.2
### ref: https://github.com/ohmeta/metapi/blob/master/metapi/Snakefile
##############################################################

import os
import sys
import pandas

shell.executable("bash")
configfile: "config.yaml"

def parse_samples(samples_tsv):
    return pandas.read_csv(samples_tsv, sep='\t').set_index("id", drop=False)

def get_sample_id(sample_df, wildcards, col):
    return sample_df.loc[wildcards.sample, [col]].dropna()[0]

_samples = parse_samples("samples/sample1.txt")
_sampleid = _samples["id"].values.tolist()

rule all:
    input:
        expand("{profile}/{sampleid}.mp4.profile", profile = config["assay"]["profile"], sampleid = _sampleid),
        expand("{profile}/{sampleid}.mp4.unclassified.profile", profile = config["assay"]["profile"], sampleid = _sampleid),
        expand("{strain}/{sampleid}/{sampleid}.pkl", strain = config["assay"]["strain"], sampleid = _sampleid),
        expand("{rmhostdir}/{sampleid}.rmhost.1.fq.gz",sampleid = _sampleid, rmhostdir=config["assay"]["rmhost"]),
        expand("{rmhostdir}/{sampleid}.rmhost.2.fq.gz",sampleid = _sampleid, rmhostdir=config["assay"]["rmhost"])



### step1 trimming & remove host reads
### To reduce disk storage usage, merge trimming and remove host together.

rule filter:
    input:
        r1 = lambda wildcards: get_sample_id(_samples, wildcards, "fq1"),
        r2 = lambda wildcards: get_sample_id(_samples, wildcards, "fq2")
    output:
        trim_r1 = temp(os.path.join(config["assay"]["trimming"], "{sample}.trimmed.1.fq.gz")),
        trim_r2 = temp(os.path.join(config["assay"]["trimming"], "{sample}.trimmed.2.fq.gz")),
        html = os.path.join(config["assay"]["trimming"], "{sample}.fastp.html"),
        json = os.path.join(config["assay"]["trimming"], "{sample}.fastp.json"),
        rmhost_r1 = protected(os.path.join(config["assay"]["rmhost"], "{sample}.rmhost.1.fq.gz")),
        rmhost_r2 = protected(os.path.join(config["assay"]["rmhost"], "{sample}.rmhost.2.fq.gz"))
    params:
        min_len = config["params"]["fastp"]["min_len"],
        ad_r1 = config["params"]["fastp"]["adapter_r1"],
        ad_r2 = config["params"]["fastp"]["adapter_r2"],
        index = config["params"]["rmhost"]["bowtie2_index"],
        env = config["envpath"]["metaphlan4"]
    threads:
        config["params"]["rmhost"]["threads"]
    log:
        fastp_log = os.path.join(config["logs"]["trimming"], "{sample}.fastp.log"),
        bowtie2_log = os.path.join(config["logs"]["rmhost"], "{sample}.rmhost.log")
    run:
        shell(
        '''
        export PATH={params.env}:$PATH

        fastp -i {input.r1} -I {input.r2} -o {output.trim_r1} -O {output.trim_r2} -w {threads} --length_required {params.min_len} --adapter_sequence={params.ad_r1} --adapter_sequence_r2={params.ad_r2} -j {output.json} -h {output.html} 2> {log.fastp_log}

        bowtie2 --end-to-end --very-sensitive -p {threads} -x {params.index} -1 {output.trim_r1} -2 {output.trim_r2} 2> {log.bowtie2_log} | samtools fastq -N -c 5 -f 12 --threads {threads} -1 {output.rmhost_r1} -2 {output.rmhost_r2} -
        ''')

### step2 filter_summary
rule seqkit_stat:
    input:
        expand("{rmhost_log_dir}/{{sample}}.rmhost.{reads}.fq.gz", rmhost_log_dir = config["assay"]["rmhost"], reads = ["1","2"])
    output:
        os.path.join(config["logs"]["rmhost"], "{sample}.rmhost.reads.summary")
    params:
        env = config["envpath"]["metaphlan4"]
    shell:
        '''
        export PATH={params.env}:$PATH
        seqkit stat {input} > {output}
        '''

rule filter_summary:
    input:
        trim = expand("{trim_res}/{sample}.fastp.json", trim_res = config["assay"]["trimming"], sample = _samples.index),
        rmhost = expand("{rmhost_res}/{sample}.rmhost.reads.summary", rmhost_res = config["logs"]["rmhost"], sample = _samples.index)
    output:
        protected(os.path.join(config["results"], "filter_summary.txt"))
    params:
        trim_summary = temp(os.path.join(config["results"], "trim_summary.txt")),
        rmhost_summary = temp(os.path.join(config["results"], "rmhost_summary.txt")),
        env = config["envpath"]["metaphlan4"]
    run:
        shell(
        '''
        export PATH={params.env}:$PATH
        python rules/filter_summary.py -t {input.trim} > {params.trim_summary}
        python rules/filter_summary.py -r {input.rmhost} > {params.rmhost_summary}
        python rules/merge_summary.py {params.trim_summary} {params.rmhost_summary} {output}
        rm {params.trim_summary} {params.rmhost_summary}
        ''')

### step3 metaphlan4
rule metaphlan4:
    input:
        r1 = os.path.join(config["assay"]["rmhost"], "{sample}.rmhost.1.fq.gz"),
        r2 = os.path.join(config["assay"]["rmhost"], "{sample}.rmhost.2.fq.gz")
    output:
        mpa = protected(os.path.join(config["assay"]["profile"], "{sample}.mp4.profile")),
        unclassified = protected(os.path.join(config["assay"]["profile"], "{sample}.mp4.unclassified.profile")),
        sam = temp(os.path.join(config["assay"]["sam"], "{sample}.sam.bz2")),
        strain = protected(os.path.join(config["assay"]["strain"], "{sample}", "{sample}.pkl"))
    threads:
        config["params"]["metaphlan4"]["threads"]
    params:
        bowtie2db = config["params"]["metaphlan4"]["bowtie2db"],
        index = config["params"]["metaphlan4"]["index"],
        bw2 = protected(os.path.join(config["assay"]["profile"], "{sample}.mp4.bw2.bz2")),
        env = config["envpath"]["metaphlan4"],
        strainpath = os.path.join(config["assay"]["strain"], "{sample}") + "/"
    shell:
        '''
        export PATH={params.env}:$PATH
        metaphlan {input.r1},{input.r2} --bowtie2db {params.bowtie2db} --index {params.index} --nproc {threads} --input_type fastq -s {output.sam} --bowtie2out {params.bw2} -t rel_ab_w_read_stats > {output.mpa}
        {params.env}/sample2markers.py -i {output.sam}  -d {params.bowtie2db}  -o {params.strainpath} -n {threads} 
        metaphlan {params.bw2} --bowtie2db {params.bowtie2db} --index {params.index} --nproc {threads} --input_type bowtie2out --unclassified_estimation > {output.unclassified}
        '''
