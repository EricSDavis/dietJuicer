# dietJuicer

Like `juicer`, but without all those calories

`dietJuicer` is a lighter-weight, HPC flexible version of juicer written with snakemake.


## Quickstart
-----------------------------------

`dietJuicer` is designed to be run on a single Hi-C library. Sequencing replicates (i.e. the same Hi-C library split to maximize complexity), can combined and run together. Run the following processing steps for each Hi-C library:

1. Clone workflow into working directory:

    ```bash
    git clone https://github.com/EricSDavis/dietJuicer.git .
    ```

2. Edit the tab-separated `samples.txt` file with paths to `Read1` and `Read2` gzipped fastq files. Multiple sequencing replicates can be combined at this step by adding a new line of paired-end files for each sequencing replicate. No naming convention is needed for fastq files, as long as they are gzipped.

    ```
    Read1   Read2
    /path/to/SeqRep1_R1.fastq.gz    /path/to/SeqRep1_R2.fastq.gz   
    /path/to/SeqRep2_R1.fastq.gz    /path/to/SeqRep2_R2.fastq.gz
    ```

3. Edit `config/config.yaml` file for your system/experiment:

    ```yaml
    ## Location of Sample-sheet containing read paths
    samples: 
        'samples.txt'

    ## Prefix to add to each output file
    prefix:
        'CI_THP1_A_1.1.1'

    ## Genome-specific reference parameters
    fasta:
        '/proj/seq/data/HG19_UCSC/Sequence/BWAIndex/genome.fa' 

    site:
        "MboI" #or 'none'

    site_file:
        'restriction_sites/hg19_MboI.txt'

    ## Rule-specific parameters
    split:
        splitsize: 2000000

    countLigations:
        ligation: "GATCGATC"

    chimera:
        mapq0_reads_included: 0
    ```

4. Submit to `SLURM` with sbatch script:

    ```bash
    sbatch runJuicer.sh
    ```

    `runJuicer.sh` will submit a long-running, low-resource job script that will spawn other jobs in the pipeline as dependencies are fulfilled. For systems other than slurm, edit `runJuicer.sh` to create a submission script for your HPC's job scheduler.

After running these steps the pipeline will produce an `output/{prefix}_dedup_merged_nodups.txt` file (and many more files), that can be used for the creating `.hic` maps.

It is recommended to fork this repo and make adjustments for your HPC environment for future Hi-C processing.


## Workflow
------------------------

The dietJuicer pipeline runs in two stages: 1) `splitFASTQ` and 2) `alignFASTQ`. Each stage runs as a separate snakemake workflow, first splitting reads according to the `splitsize` parameter in the `config/config.yaml` file, and then following the traditional `juicer` pipeline steps in `alignFASTQ`. The pipeline results in `stats` files, and a `merged_nodups` file which can be used to create a `.hic` file.

See the diagrams below for a DAG representation of the workflow:

### Step 1: splitFASTQ subworkflow

![](img/splitFASTQdag.png)


### Step 2: alignFASTQ subworkflow

![](img/alignFASTQdag.png)


## Setup & Dependencies

* snakemake version 5.10.0
    ```bash
    python3 -m venv env
    source minenv/bin/activate
    pip3 install --upgrade snakemake==5.10.0
    ```

* python 3 package `psutil` for benchmarking
    ```
    module load python/3.6.6
    pip3.6 install --user psutil
    ```