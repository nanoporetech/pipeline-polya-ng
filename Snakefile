
import os
from os import path

if not workflow.overwrite_configfiles:
    configfile: "config.yml"
workdir: path.join(config["workdir_top"], config["pipeline"])

WORKDIR = path.join(config["workdir_top"], config["pipeline"])
SNAKEDIR = path.dirname(workflow.snakefile)

if(config["per_transcript_plots"]):
    config["per_transcript_plots"] = "true"

include: "snakelib/utils.snake"


rule dump_versions:
    output:
        ver = "versions.txt"
    conda: "env.yml"
    shell:"""
    command -v conda > /dev/null && conda list > {output.ver}
    """

rule make_reference:
    input:
        trs = config["transcriptome"],
    output:
        ref = "input/reference.fas"
    params:
        spike = config["spikein_fasta"],
    shell:"""
    cp {input.trs} {output.ref}
    if [ ! -z "{params.spike}" ]
    then
        cat {params.spike} >> {output.ref}
    fi
    """

rule make_fofn:
    input:
        sum_dir = config["summary_dir"]
    output:
        fofn = "input/summaries.fofn"
    shell:"""
    find {input.sum_dir} -maxdepth 1 -type f -name "sequencing_summary*.txt" -print > {output.fofn}
    """

rule make_fastq:
    input:
        fqd = config["fastq_dir"]
    output:
        fastq = "input/reads.fastq"
    conda: "env.yml"
    shell:"""
    find {input.fqd} -maxdepth 1 -type f -name "*.fastq" -exec cat {{}} \; | seqkit seq --rna2dna - | seqkit rmdup -n - > {output.fastq} 
    """

rule nanopolish_index:
    input:
        f5dir = config["fast5_dir"],
        fq = rules.make_fastq.output.fastq,
        fofn = rules.make_fofn.output.fofn,
    output:
        fastq_index = "input/reads.fastq.index",
        fastq_index_fai = "input/reads.fastq.index.fai",
        fastq_index_gzi = "input/reads.fastq.index.gzi",
        fastq_index_readdb = "input/reads.fastq.index.readdb",
    conda: "env.yml"
    shell:"""
    nanopolish index -d {input.f5dir} -f {input.fofn} {input.fq}
    """

rule map_reads:
    input:
        fastq = rules.make_fastq.output.fastq,
        ref = rules.make_reference.output.ref,
    output:
        aln = "alignment/aligned_reads_sorted.bam",
        alni = "alignment/aligned_reads_sorted.bam.bai",
    params:
        mmq = config["min_mapping_qual"]
    threads: config["threads"]
    conda: "env.yml"
    shell:"""
    minimap2 -ax map-ont -t {threads} {input.ref} {input.fastq} | samtools view -b -F 2304 -q {params.mmq} | samtools sort -@ {threads} -o {output.aln};
    samtools index {output.aln}
    """

rule call_tails:
    input:
        index = rules.nanopolish_index.output.fastq_index,
        fastq = rules.make_fastq.output.fastq,
        bam = rules.map_reads.output.aln,
        ref = rules.make_reference.output.ref,
    output:
        tails = "tails/all_tails.tsv"
    threads: config["threads"]
    conda: "env.yml"
    shell:"""
    nanopolish polya -r {input.fastq} -b {input.bam} -g {input.ref} -t {threads} > {output.tails}
    """

rule filter_tails:
    input:
        tails = rules.call_tails.output.tails,
    output:
        flt = temp("tails/filtered_tails_sp.tsv"),
        pdf = "reports/filtering_report.pdf",
        tsv = "reports/filtering_report.tsv",
    conda: "env.yml"
    shell:"""
    {SNAKEDIR}/scripts/filter_tails.py -i {input.tails} -o {output.flt} -r {output.pdf} -t {output.tsv}
    """

rule spikein_qc:
    input:
        tails = rules.filter_tails.output.flt,
    output:
        tsv = "tails/spikein_tails.tsv",
        pdf = "reports/spikein_report.pdf",
        med = "reports/spikein_medians.tsv",
        dat = "tails/filtered_tails.tsv",
    params:
        spike = config["spikein_fasta"]
    conda: "env.yml"
    shell:"""
    if [ ! -z "{params.spike}" ]
    then
        {SNAKEDIR}/scripts/spikein_qc.py -i {input.tails} -o {output.tsv} -d {output.dat} -m {output.med} -r {output.pdf} -s {params.spike}
    else
        touch {output.tsv} {output.pdf} {output.med}
        mv {input.tails} {output.dat}
    fi
    """

rule tails_report:
    input:
        tails = rules.spikein_qc.output.dat,
    output:
        pdf = "reports/tails_report.pdf",
        tsv = "reports/tails_report.tsv",
    params:
        x = config["per_transcript_plots"]
    conda: "env.yml"
    shell:"""
    X=""
    if [ {params.x} == "true" ];
    then
        X="-x"
    fi
    {SNAKEDIR}/scripts/tails_report.py $X -i {input.tails} -m {output.tsv} -r {output.pdf}
    """

rule all:
    input:
        ver = rules.dump_versions.output.ver, 
        ref = rules.make_reference.output.ref,
        sum_fofn = rules.make_fofn.output.fofn,
        fastq = rules.make_fastq.output.fastq,
        index = rules.nanopolish_index.output.fastq_index,
        aln = rules.map_reads.output.aln,
        tails = rules.call_tails.output.tails,
        sp_rep = rules.spikein_qc.output.pdf,
        tails_rep = rules.tails_report.output.pdf,
