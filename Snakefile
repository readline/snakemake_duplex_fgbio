import os,sys
import yaml
from os.path import join
import pandas as pd
from scripts.load import samplesheet
from scripts.utils import (allocated, ignore)

# Load config file
configfile: 'config.yaml'
workdir: config['workdir']

pipedir = config['pipelinedir']
print(config)

# Load cluster config
with open(join(config['pipelinedir'], 'cluster.yaml')) as infile:
    cluster = yaml.safe_load(infile)
print(cluster)

# Load sample sheet
sampledic, rundic, run2sample = samplesheet(join(config['pipelinedir'],'samplesheet.tsv'))

print(sampledic)
print(rundic)
print(run2sample)


rule all:
    input: 
        expand(
            join(config['workdir'], "01.MergeData", "{sample}", "{sample}.QC.html"),
            sample=sampledic,
        ),

        expand(
            join(config['workdir'], "02.Duplex", "{sample}", "{sample}.cons.mapped.bam"),
            sample=sampledic,
        ),

        expand(
            join(config['workdir'], "02.Duplex", "{sample}", "{sample}.duplex_seq_metrics.umi_counts.txt"),
            sample=sampledic,
        ),

rule Merge1:
    input:
        reads1 = lambda wildcards: [rundic[rg]['r1'] for rg in sampledic[wildcards.sample]],
    output:
        reads1out = temp(join(config['workdir'], "01.MergeData", "{sample}", "{sample}.R1.fq.gz")),
    log: 
        out = join(config['pipelinedir'], "logs", "Merge1", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "Merge1", "{sample}.e"),
    threads:
        int(allocated("threads", "Merge1", cluster))
    shell:
        """
        if [ $(ls -1 {input.reads1} | wc -l) -eq 1 ]; then
            ln -s {input.reads1} {output.reads1out} \
               >> {log.out} 2> {log.err}
        else
            zcat {input.reads1} | gzip -c - > {output.reads1out} 2> {log.err}
        fi
        """

rule Merge2:
    input:
        reads2 = lambda wildcards: [rundic[rg]['r2'] for rg in sampledic[wildcards.sample]],
    output:
        reads2out = temp(join(config['workdir'], "01.MergeData", "{sample}", "{sample}.R2.fq.gz")),
    log: 
        out = join(config['pipelinedir'], "logs", "Merge2", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "Merge2", "{sample}.e"),
    threads:
        int(allocated("threads", "Merge2", cluster))
    shell:
        """
        if [ $(ls -1 {input.reads2} | wc -l) -eq 1 ]; then
            ln -s {input.reads2} {output.reads2out} \
               >> {log.out} 2> {log.err}
        else
            zcat {input.reads2} | gzip -c - > {output.reads2out} 2> {log.err}
        fi
        """

rule Fastp:
    input:
        reads1 = join(config['workdir'], "01.MergeData", "{sample}", "{sample}.R1.fq.gz"),
        reads2 = join(config['workdir'], "01.MergeData", "{sample}", "{sample}.R2.fq.gz"),
    output:
        reads1 = join(config['workdir'], "01.MergeData", "{sample}", "{sample}.R1.cln.fq.gz"),
        reads2 = join(config['workdir'], "01.MergeData", "{sample}", "{sample}.R2.cln.fq.gz"),
        htmlout = join(config['workdir'], "01.MergeData", "{sample}", "{sample}.QC.html"),
        jsonout = join(config['workdir'], "01.MergeData", "{sample}", "{sample}.QC.json"),
    log:
        out = join(config['pipelinedir'], "logs", "Fastp", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "Fastp", "{sample}.e"),
    threads:
        int(allocated("threads", "Fastp", cluster))
    container:
        config['container']['fastp']
    shell:
        '''
        fastp -i {input.reads1} \
            -I {input.reads2} \
            -o {output.reads1} \
            -O {output.reads2} \
            -h {output.htmlout} \
            -j {output.jsonout} \
            -l 30 \
            -w {threads} \
            >> /dev/null 2> {log.err}
        '''

rule FastqToBam:
    input:
        reads1 = join(config['workdir'], "01.MergeData", "{sample}", "{sample}.R1.cln.fq.gz"),
        reads2 = join(config['workdir'], "01.MergeData", "{sample}", "{sample}.R2.cln.fq.gz"),
    output:
        bam    = join(config['workdir'], "02.Duplex", "{sample}", "{sample}.unmapped.bam"),
    params:
        structure = config['structure']
    log:
        out = join(config['pipelinedir'], "logs", "FastqToBam", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "FastqToBam", "{sample}.e"),
    threads:
        int(allocated("threads", "FastqToBam", cluster))
    container:
        config['container']['fgbio']
    shell:
        '''
        fgbio \
          -Xmx6g \
          --tmp-dir=. \
          --async-io=true \
          --compression=1 \
          FastqToBam \
            --input {input.reads1} {input.reads2} \
            --output {output.bam} \
            --read-structures {params.structure} \
            --sample {wildcards.sample} \
            --library {wildcards.sample}
        '''

rule Align_unmapped:
    input:
        bam    = join(config['workdir'], "02.Duplex", "{sample}", "{sample}.unmapped.bam"),
    output:
        bam    = join(config['workdir'], "02.Duplex", "{sample}", "{sample}.mapped.bam"),
    log:
        out = join(config['pipelinedir'], "logs", "04.Align_unmapped", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "04.Align_unmapped", "{sample}.e"),
    threads:
        int(allocated("threads", "Align_unmapped", cluster))
    container:
        config['container']['fgbio']
    shell:
        '''
        samtools \
            fastq {input.bam} \
            2>> {log.err} \
            | \
        bwa mem \
            -t {threads} \
            -p \
            -K 150000000 \
            -Y {config[cachedir]}/{config[references][gatkbundle][path]}/Homo_sapiens_assembly38.fasta \
            - \
            2>> {log.err} \
            | \
        fgbio \
            -Xmx8g \
            --async-io=true \
            --compression=0 \
            ZipperBams \
            --unmapped {input.bam} \
            --ref {config[cachedir]}/{config[references][gatkbundle][path]}/Homo_sapiens_assembly38.fasta \
            --output /dev/stdout \
            --tags-to-reverse Consensus \
            --tags-to-revcomp Consensus \
            2>> {log.err} | \
        samtools \
            sort \
            --template-coordinate \
            --threads 4 \
            -o {output.bam} \
            >> {log.out} 2>> {log.err}
        '''

rule GroupReadsByUmi:
    input:
        bam    = join(config['workdir'], "02.Duplex", "{sample}", "{sample}.mapped.bam"),
    output:
        bam    = join(config['workdir'], "02.Duplex", "{sample}", "{sample}.grouped.bam"),
        hist   = join(config['workdir'], "02.Duplex", "{sample}", "{sample}.grouped-family-sizes.txt")
    log:
        out = join(config['pipelinedir'], "logs", "05.GroupReadsByUmi", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "05.GroupReadsByUmi", "{sample}.e"),
    threads:
        int(allocated("threads", "GroupReadsByUmi", cluster))
    container:
        config['container']['fgbio']
    shell:
        '''
        fgbio \
            -Xmx8g \
            --async-io=true \
            --compression=1 \
            GroupReadsByUmi \
            --strategy Paired \
            --edits 1 \
            --input {input.bam} \
            --output {output.bam} \
            --family-size-histogram {output.hist} \
            --min-map-q=10 \
            >> {log.out} 2>> {log.err}
        '''

rule CallDuplexConsensusReads:
    input:
        bam    = join(config['workdir'], "02.Duplex", "{sample}", "{sample}.grouped.bam"),
    output:
        bam    = join(config['workdir'], "02.Duplex", "{sample}", "{sample}.cons.unmapped.bam"),
    log:
        out = join(config['pipelinedir'], "logs", "06.CallDuplexConsensusReads", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "06.CallDuplexConsensusReads", "{sample}.e"),
    threads:
        int(allocated("threads", "CallDuplexConsensusReads", cluster))
    container:
        config['container']['fgbio']
    shell:
        '''
        fgbio \
            -Xmx8g \
            --async-io=true \
            --compression=1 \
            CallDuplexConsensusReads \
            --input {input.bam} \
            --output {output.bam} \
            --min-reads 2 1 1 \
            --min-input-base-quality 30 \
            --threads {threads} \
            >> {log.out} 2>> {log.err}
        '''

rule CollectDuplexSeqMetrics:
    input:
        bam    = join(config['workdir'], "02.Duplex", "{sample}", "{sample}.grouped.bam"),
    output:
        txt    = join(config['workdir'], "02.Duplex", "{sample}", "{sample}.duplex_seq_metrics.umi_counts.txt"),
    params:
        prefix = join(config['workdir'], "02.Duplex", "{sample}", "{sample}.duplex_seq_metrics"),
    log:
        out = join(config['pipelinedir'], "logs", "07.CollectDuplexSeqMetrics", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "07.CollectDuplexSeqMetrics", "{sample}.e"),
    threads:
        int(allocated("threads", "CollectDuplexSeqMetrics", cluster))
    container:
        config['container']['fgbio']
    shell:
        '''
        fgbio \
            -Xmx8g \
            --async-io=true \
            --compression=1 \
            CollectDuplexSeqMetrics \
            --input {input.bam} \
            --output {params.prefix} \
            --duplex-umi-counts=true \
            >> {log.out} 2>> {log.err}
        '''

rule Align_consensus:
    input:
        bam    = join(config['workdir'], "02.Duplex", "{sample}", "{sample}.cons.unmapped.bam"),
    output:
        bam    = join(config['workdir'], "02.Duplex", "{sample}", "{sample}.cons.mapped.bam"),
    log:
        out = join(config['pipelinedir'], "logs", "08.Align_consensus", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "08.Align_consensus", "{sample}.e"),
    threads:
        int(allocated("threads", "Align_consensus", cluster))
    container:
        config['container']['fgbio']
    shell:
        '''
        samtools \
            fastq \
            {input.bam} \
            2>> {log.err} \
            | \
        bwa \
            mem \
            -t {threads} \
            -p \
            -K 150000000 \
            -Y {config[cachedir]}/{config[references][gatkbundle][path]}/Homo_sapiens_assembly38.fasta \
            - \
            2>> {log.err} \
            | \
        fgbio \
            -Xmx8g \
            --async-io=true \
            --compression=1 \
        ZipperBams \
            --unmapped {input.bam} \
            --ref {config[cachedir]}/{config[references][gatkbundle][path]}/Homo_sapiens_assembly38.fasta \
            --output {output.bam} \
            --tags-to-reverse Consensus \
            --tags-to-revcomp Consensus \
            >> {log.out} 2>> {log.err}
        '''
