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
		expand('output/{prefix}_countLigations_split{sample}_norm.txt.res.txt', prefix=config["prefix"], sample=samples.index),
		expand('output/{prefix}_countLigations_split{sample}_linecount.txt', prefix=config["prefix"], sample=samples.index)

rule countLigations:
	input:
		R1 = lambda wildcards: samples.iloc[[int(i) for i in wildcards.sample]]["Read1"],
		R2 = lambda wildcards: samples.iloc[[int(i) for i in wildcards.sample]]["Read2"]
	output:
		temp = temp('output/{prefix}_countLigations_split{sample}_temp'),
		res = 'output/{prefix}_countLigations_split{sample}_norm.txt.res.txt',
		linecount = 'output/{prefix}_countLigations_split{sample}_linecount.txt'
	params:
		ligation = config['countLigations']['ligation']
	run:
		shell("R1={input.R1} R2={input.R2} ligation={params.ligation} temp={output.temp} res={output.res} linecount={output.linecount} ./scripts/countLigations.sh")