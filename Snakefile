#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import glob

##### Load config and sample sheets #####

configfile: "config/config.yaml"

samples = pd.read_table(config["samples"])

##### Error handlers ####

onsuccess:
	print("Workflow completed successfully!")

##### Define rules #####
rule all:
	input:
		expand('output/{prefix}makeFile_split{sample}_R{pair}.txt', sample=list(samples.index), pair=[1,2], prefix=config["prefix"])

rule makeFiles:
	input:
		in1 = lambda wildcards: samples["Read1"],
		in2 = lambda wildcards: samples["Read2"]
	output:
		out1 = 'output/{prefix}makeFile_split{sample}_R1.txt',
		out2 = 'output/{prefix}makeFile_split{sample}_R2.txt'
	message: "Launching makeFiles job"
	shadow: "full"
	shell:
		'touch {output.out1} {output.out2};'