#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import glob

configfile: "config/config.yaml"

## Success Message for pipeline
onsuccess:
	print("Entire workflow completed successfully!")

## Split samples
subworkflow splitFASTQ:
	workdir: '.'
	snakefile: 'splitFASTQ'
	configfile: 'config/config.yaml'

## Align splits and create merged_nodups
subworkflow alignFASTQ:
	workdir: '.'
	snakefile: 'alignFASTQ'
	configfile: 'config/config.yaml'

## Trigger each subworkflow
rule all:
	input:
		splitFASTQ('output/splitFASTQ_done.txt'),
		alignFASTQ(expand('output/{prefix}_merge_{extension}.bam', prefix=config["prefix"], extension=['collisions', 'collisions_low_mapq', 'unmapped', 'mapq0'])),
		alignFASTQ(expand('output/{prefix}_dedupAlignablePart3_alignable.bam', prefix=config["prefix"])),
		alignFASTQ(expand('output/{prefix}_dupStats_stats_dups{extension}', prefix=config["prefix"], extension=['.txt', '_hists.m'])),
		alignFASTQ(expand('output/{prefix}_inter.txt', prefix=config["prefix"])),
		alignFASTQ(expand('output/{prefix}_inter_30.txt', prefix=config["prefix"]))