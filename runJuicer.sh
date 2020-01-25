#! /bin/bash -login
#SBATCH -J scheduler
#SBATCH -t 7-00:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1  
#SBATCH -p general
#SBATCH --mem=2gb

## Exit if any command fails
set -e

## Load required modules
module load python/3.6.6
module load bwa
module load samtools

## Create and activate virtual environment with requirements
python3 -m venv env && source env/bin/activate && pip3 install -r requirements.txt

## Make directory for slurm logs
mkdir -p logs_slurm

## Execute splitFASTQ snakemake workflow
snakemake -j 50 -s splitFASTQ --latency-wait 500 --cluster-config "config/cluster.yaml" --cluster "sbatch --job-name {cluster.name} -p {cluster.partition} -n {cluster.tasks} -N {cluster.nodes} --mem={cluster.mem} -t {cluster.time} --output {cluster.output} --error {cluster.error}"

## Execute alignFASTQ snakemake workflow
snakemake -j 50 -s alignFASTQ --latency-wait 500 --cluster-config "config/cluster.yaml" --cluster "sbatch --job-name {cluster.name} -p {cluster.partition} -n {cluster.tasks} -N {cluster.nodes} --mem={cluster.mem} -t {cluster.time} --output {cluster.output} --error {cluster.error}"

## Success message
echo "Entire workflow completed successfully!"