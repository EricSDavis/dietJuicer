#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import glob
import shutil

##### Load config and sample sheets #####
configfile: "config/config.yaml"

## Read in samplesheet
samples = pd.read_table(config["samplesheet"])

## Convert all columns to strings
samples = samples.astype(str)

## Concatenate the sequencing directory to Read1 and Read2 for full paths
samples['Read1'] = samples[['Sequencing_Directory', 'Read1']].apply(lambda row: os.path.join(*row), axis=1)
samples['Read2'] = samples[['Sequencing_Directory', 'Read2']].apply(lambda row: os.path.join(*row), axis=1)

## Concatenate columns to identify which groups to run (i.e. Seq_Rep will be run together)
samples['id'] = samples[config['mergeBy']].agg('_'.join, axis=1)

## Group by id and extract Read1 & Read2
read1 = samples.groupby('id')['Read1'].apply(list).to_dict()
read2 = samples.groupby('id')['Read2'].apply(list).to_dict()

## Get wildcards for splitNames
## Build {group:splitName} dictionary for read1 from completed split file (i.e. 'output/{group}/splitR1_done.txt')
splitRead1 = dict()
for i in read1.keys():
    filename1 = os.path.join('output', i, '{}_splitR1_done.txt'.format(i))
    with open(filename1) as f:
        split1 = f.read().splitlines()
        splitNames1 = [j[0:4] for j in split1]
    splitRead1[i] = splitNames1

## Build {group:splitName} dictionary for read2 from completed split file (i.e. 'output/{group}/splitR2_done.txt')
splitRead2 = dict()
for i in read2.keys():
    filename2 = os.path.join('output', i, '{}_splitR2_done.txt'.format(i))
    with open(filename2) as f:
        split2 = f.read().splitlines()
        splitNames2 = [j[0:4] for j in split2]
    splitRead2[i] = splitNames2

## Build splitNames = {group:splitName} dictionary if both dictionaries match
if splitRead1 == splitRead2:
    splitNames = splitRead1
else:
    splitNames = dict()

## Define actions on success
onsuccess:
    ## Success message
    print("alignFASTQ completed successfully!")

    ## Remove split files and folders
    for key in splitNames:
        shutil.rmtree(os.path.join('output', key, 'splitsR1'))
        shutil.rmtree(os.path.join('output', key, 'splitsR2'))
        os.remove(os.path.join('output', key, '{group}_splitR1_done.txt'.format(group=key)))
        os.remove(os.path.join('output', key, '{group}_splitR2_done.txt'.format(group=key)))
        print("output/{group}/splitsR1 and output/{group}/splitsR2 have been removed".format(group=key))


##### Define rules #####
rule all:
    input:
        [expand("output/{group}/{group}_dedup_merged_nodups.txt.gz", group=key) for key in splitNames]
      
rule countLigations:
    input:
        R1 = lambda wildcards: ['output/{group}/splitsR1/{splitName}_R1.fastq.gz'.format(group=wildcards.group, splitName=wildcards.splitName)],
        R2 = lambda wildcards: ['output/{group}/splitsR2/{splitName}_R2.fastq.gz'.format(group=wildcards.group, splitName=wildcards.splitName)]
    output:
        temp = temp('output/{group}/{group}_countLigations_split{splitName}_temp'),
        res = temp('output/{group}/{group}_countLigations_split{splitName}_norm.txt.res.txt'),
        linecount = temp('output/{group}/{group}_countLigations_split{splitName}_linecount.txt')
    log:
        err = "output/{group}/logs/{group}_countLigations_split{splitName}.err",
        out = "output/{group}/logs/{group}_countLigations_split{splitName}.out"
    params:
        ligation = config['ligation']
    threads: 1
    benchmark: 
        'output/{group}/benchmarks/{group}_countLigations_split{splitName}.tsv'
    shell:
        "R1={input.R1} R2={input.R2} ligation={params.ligation} temp={output.temp} res={output.res} linecount={output.linecount} sh scripts/countLigations.sh 2> {log.err} 1> {log.out}"

rule align:
    input:
        R1 = lambda wildcards: ['output/{group}/splitsR1/{splitName}_R1.fastq.gz'.format(group=wildcards.group, splitName=wildcards.splitName)],
        R2 = lambda wildcards: ['output/{group}/splitsR2/{splitName}_R2.fastq.gz'.format(group=wildcards.group, splitName=wildcards.splitName)]
    output:
        sam = temp("output/{group}/{group}_align_split{splitName}.sam")
    log:
        err = "output/{group}/logs/{group}_align_split{splitName}.err"
    threads: 8
    params:
        fasta = config['fasta']
    benchmark: 
        "output/{group}/benchmarks/{group}_align_split{splitName}.tsv"
    shell:
        "module load bwa; "
        "bwa mem -SP5M -t {threads} {params.fasta} {input.R1} {input.R2} > {output.sam} 2> {log.err}"

rule chimera:
    input:
        sam = rules.align.output.sam
    output:
        norm = temp("output/{group}/{group}_align_split{splitName}_norm.txt"),
        normRes = temp("output/{group}/{group}_align_split{splitName}_norm.txt.res.txt"),
        alignable = temp("output/{group}/{group}_align_split{splitName}_alignable.sam"),
        collisions = temp("output/{group}/{group}_align_split{splitName}_collisions.sam"),
        lowqcollisions = temp("output/{group}/{group}_align_split{splitName}_collisions_low_mapq.sam"),
        unmapped = temp("output/{group}/{group}_align_split{splitName}_unmapped.sam"),
        mapq0 = temp("output/{group}/{group}_align_split{splitName}_mapq0.sam"),
        unpaired = temp("output/{group}/{group}_align_split{splitName}_unpaired.sam"),
        done = temp("output/{group}/{group}_chimera_split{splitName}_done.txt")
    log:
        err = "output/{group}/logs/{group}_chimera_split{splitName}.err",
        out = "output/{group}/logs/{group}_chimera_split{splitName}.out"
    threads: 1
    params:
        fname = 'output/{group}/{group}_align_split{splitName}',
        mapq0_reads_included = config['mapq0_reads_included']
    benchmark: 
        'output/{group}/benchmarks/{group}_chimera_split{splitName}.tsv'
    shell:
        """
        touch {output.norm} {output.normRes} {output.alignable} {output.collisions} {output.lowqcollisions} {output.unmapped} {output.mapq0} {output.unpaired}
        gawk -v "fname"={params.fname} -v "mapq0_reads_included"={params.mapq0_reads_included} -f ./scripts/chimeric_blacklist.awk {input.sam} 2> {log.err} 1> {log.out} && \
        touch {output.done}
        """

rule fragment:
    input:
        done = rules.chimera.output.done,
        norm = rules.chimera.output.norm
    output:
        frag = temp("output/{group}/{group}_fragment_split{splitName}.frag.txt")
    log:
        err = "output/{group}/logs/{group}_fragment_split{splitName}.err"
    threads: 1
    params:
        site = config['site'],
        site_file = config['site_file']
    benchmark:
        'output/{group}/benchmarks/{group}_fragment_split{splitName}.tsv'
    shell:
        """
        if [ {params.site} != "none" ]
        then
                ./scripts/fragment.pl {input.norm} {output.frag} {params.site_file} 2> {log.err}    
        else
                awk '{{printf("%s %s %s %d %s %s %s %d", $1, $2, $3, 0, $4, $5, $6, 1); for (i=7; i<=NF; i++) {{printf(" %s",$i);}}printf("\\n");}}' {input.norm} > {output.frag} 2> {log.err}
        fi
        """

rule sort:
    input:
        rules.fragment.output.frag
    output:
        sorted = temp("output/{group}/{group}_sort_split{splitName}.sort.txt")
    log:
        err = "output/{group}/logs/{group}_sort_split{splitName}.err"
    threads: 8
    shadow: "minimal"
    benchmark:
        'output/{group}/benchmarks/{group}_sort_split{splitName}.tsv'
    shell:
        """
        tmp_dir=$(mktemp -d -p $PWD)
        sort --parallel={threads} -T $tmp_dir -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n {input} > {output.sorted} 2> {log.err}
        rm -rf $tmp_dir
        """

rule mergedSort:
    input:
        lambda wildcards: ['output/{group}/{group}_sort_split{splitName}.sort.txt'.format(group=wildcards.group, splitName=value) for value in splitNames[wildcards.group]]
    output:
        temp('output/{group}/{group}_mergedSort_merged_sort.txt')
    log:
        err = 'output/{group}/logs/{group}_mergedSort.err'
    threads: 8
    shadow: "minimal"
    benchmark:
        'output/{group}/benchmarks/{group}_mergedSort.tsv'
    shell:
        """
        tmp_dir=$(mktemp -d -p $PWD)
        sort --parallel={threads} -T $tmp_dir -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n {input} > {output} 2> {log.err}
        rm -rf $tmp_dir
        """

rule dedup:
    input:
        rules.mergedSort.output
    output:
        dups = temp("output/{group}/{group}_dedup_dups.txt"),
        merged_nodups = temp("output/{group}/{group}_dedup_merged_nodups.txt"),
        optdups = temp("output/{group}/{group}_dedup_opt_dups.txt"),
        done = temp('output/{group}/{group}_dedup_done.txt')
    log:
        err = 'output/{group}/logs/{group}_dedup.err'
    params:
        name = 'output/{group}/{group}_'
    threads: 1
    benchmark:
        'output/{group}/benchmarks/{group}_dedup.tsv'
    shell:
        """
        touch {output.dups} {output.merged_nodups} {output.optdups}
        awk -f ./scripts/dups.awk -v name={params.name} {input} 2> {log.err} && \
        touch {output.done}
        """

def interInput(wildcards):
    combined=[]
    for i in splitNames[wildcards.group]:
        combined.append('output/{group}/{group}_countLigations_split{splitName}_norm.txt.res.txt'.format(group=wildcards.group, splitName=i))
        combined.append('output/{group}/{group}_align_split{splitName}_norm.txt.res.txt'.format(group=wildcards.group, splitName=i))
    return combined

rule inter:
    input:
        done = rules.dedup.output.done,
        res = interInput,
        dups = rules.dedup.output.dups,
        optdups = rules.dedup.output.optdups,
        merged_nodups = rules.dedup.output.merged_nodups
    output:
        inter = 'output/{group}/{group}_inter.txt',
        hists = 'output/{group}/{group}_inter_hists.m'
    log:
        err = 'output/{group}/logs/{group}_inter.err'
    params:
        ligation = config['ligation'],
        site_file = config['site_file']
    threads: 1
    benchmark:
        'output/{group}/benchmarks/{group}_inter.tsv'
    shell:
        """
        cat {input.res} | awk -f ./scripts/stats_sub.awk > {output.inter} 2> {log.err}
        ./scripts/juicer_tools LibraryComplexity $(wc -l < {input.merged_nodups}) $(wc -l < {input.dups}) $(wc -l < {input.optdups}) >> {output.inter} 2>> {log.err}
        ./scripts/statistics.pl -s {params.site_file} -l {params.ligation} -o {output.inter} -q 1 {input.merged_nodups} 2>> {log.err}
        """

rule inter30:
    input:
        done = rules.dedup.output.done,
        res = interInput,
        dups = rules.dedup.output.dups,
        optdups = rules.dedup.output.optdups,
        merged_nodups = rules.dedup.output.merged_nodups
    output:
        inter30 = 'output/{group}/{group}_inter_30.txt',
        hists30 = 'output/{group}/{group}_inter_30_hists.m'
    log:
        err = 'output/{group}/logs/{group}_inter30.err'
    params:
        ligation = config['ligation'],
        site_file = config['site_file']
    threads: 1
    benchmark:
        'output/{group}/benchmarks/{group}_inter30.tsv'
    shell:
        """
        cat {input.res} | awk -f ./scripts/stats_sub.awk > {output.inter30} 2>> {log.err}
        ./scripts/juicer_tools LibraryComplexity $(wc -l < {input.merged_nodups}) $(wc -l < {input.dups}) $(wc -l < {input.optdups}) >> {output.inter30} 2>> {log.err}
        ./scripts/statistics.pl -s {params.site_file} -l {params.ligation} -o {output.inter30} -q 30 {input.merged_nodups} 2>> {log.err}
        """

rule hic:
    input:
        done = rules.dedup.output.done,
        inter = rules.inter.output.inter,
        hists = rules.inter.output.hists,
        merged_nodups = rules.dedup.output.merged_nodups
    output:
        hic = 'output/{group}/{group}_inter.hic'
    log:
        err = 'output/{group}/logs/{group}_hic.err',
        out = 'output/{group}/logs/{group}_hic.out'
    params:
        chromSizes = config['chromSizes']
    threads: 1
    benchmark:
        'output/{group}/benchmarks/{group}_hic.tsv'
    shell:
        """
        ./scripts/juicer_tools pre -s {input.inter} -g {input.hists} -q 1 {input.merged_nodups} {output.hic} {params.chromSizes} 1> {log.out} 2> {log.err}
        """

rule hic30:
    input:
        done = rules.dedup.output.done,
        inter30 = rules.inter30.output.inter30,
        hists30 = rules.inter30.output.hists30,
        merged_nodups = rules.dedup.output.merged_nodups
    output:
        hic30 = 'output/{group}/{group}_inter_30.hic'
    log:
        err = 'output/{group}/logs/{group}_hic30.err',
        out = 'output/{group}/logs/{group}_hic30.out'
    params:
        chromSizes = config['chromSizes']
    threads: 1
    benchmark:
        'output/{group}/benchmarks/{group}_hic30.tsv'
    shell:
        """
        ./scripts/juicer_tools pre -s {input.inter30} -g {input.hists30} -q 30 {input.merged_nodups} {output.hic30} {params.chromSizes} 1> {log.out} 2> {log.err}
        """

rule compress:
    input:
        done = rules.dedup.output.done,
        hic = rules.hic.output.hic, 
        hic30 = rules.hic30.output.hic30,
        merged_nodups = rules.dedup.output.merged_nodups
    output:
        'output/{group}/{group}_dedup_merged_nodups.txt.gz'
    log:
        err = 'output/{group}/logs/{group}_compress.err'
    threads: 8
    benchmark:
        'output/{group}/benchmarks/{group}_compress.tsv'
    shell:
        """
        module load pigz
        pigz -p {threads} {input.merged_nodups}
        """