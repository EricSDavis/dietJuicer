#! /bin/bash -login
#SBATCH -J scheduler
#SBATCH -t 7-00:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1  
#SBATCH -p general
#SBATCH --mem=2gb

## Load required modules
module load python/3.6.6
module load bwa
module load samtools

## Navigate to snakemake pipeline
# cd /path/to/snakemake/pipeline

## Make directory for slurm logs
mkdir -p logs_slurm

## Execute snakemake
# snakemake -s Snakefile --cluster-config cluster.yaml --cluster 'sbatch -t {cluster.time} --mem={cluster.mem} -c {cluster.tasks} -o {cluster.output} -e {cluster.error}'
snakemake -j 50 -s Snakefile --configfile config/config.yaml --latency-wait 500 --cluster-config config/cluster.yaml --cluster "sbatch --job-name {cluster.name} -p {cluster.partition} -n {cluster.tasks} -N {cluster.nodes} --mem={cluster.mem} -t {cluster.time} --output {cluster.output} --error {cluster.error}"
