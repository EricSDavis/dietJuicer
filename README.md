# dietJuicer

Like `juicer`, but without all those calories

`dietJuicer` is a lighter-weight, HPC flexible version of juicer written with snakemake.


## Quickstart
-----------------------------------

1. Clone workflow into working directory:

    For running locally use -b dietJuicerLocal
    ```bash
    git clone -b dietJuicerLocal https://github.com/EricSDavis/dietJuicer.git
    ```   

    For running on an HPC use -b dietJuicerCluster
    ```bash
    git clone -b dietJuicerCluster https://github.com/EricSDavis/dietJuicer.git
    ```

2. Edit `samples.txt` to reflect your sample information

3. Edit `config/config.yaml` file for your system/experiment:

    ```yaml
    ## Location of Sample-sheet and read paths
    samples: 
        'samples.txt'

    prefix:
        'CI_THP1_A_1.1.1'

    ## Genome-specific reference parameters
    fasta:
        '/Users/eric/Phanstiel\ Lab\ Dropbox/Eric\ Davis/genomes/hg19/bwa/hg19.fa' #'/proj/seq/data/HG19_UCSC/Sequence/BWAIndex/genome.fa' 

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

4. Launch snakemake pipeline directly or with the `runJuicer.sh` script. If running locally, `sh runJuicer.sh` will launch snakemake directly. If running on a cluster `sbatch runJuicer.sh` will submit a long-running, low resource job that will spawn other jobs in the pipeline using SLURM. For other job schedulers, copy the snakemake command and wrap it in the appropriate script.

## Setup & Dependencies

* python 3 package `psutil` for benchmarking
    ```
    module load python/3.6.6
    pip3.6 install --user psutil
    ```


## Workflow
------------------------

The dietJuicer pipeline runs in two stages: 1) `splitFASTQ` and 2) `alignFASTQ`. Each stage runs as a subworkflow, first splitting reads according to the `splitsize` parameter in the `config/config.yaml` file, and then following the traditional `juicer` pipeline steps in `alignFASTQ`. The pipeline results in `stats` files, and a `merged_nodups` file which can be used to create a `.hic` file.

See the diagrams below for a DAG representation of the workflow:

### Step 1: splitFASTQ subworkflow

![](img/splitFASTQdag.png)


### Step 2: alignFASTQ subworkflow

![](img/alignFASTQdag.png)