#!/bin/bash

## Confirm that dietJuicer Error jobs have finished correctly ####
## Slurm log to search
slurmLog="$1"

## Get line number for the end of splitFASTQ
endSplit=$(grep -n 'splitFASTQ completed successfully!' $slurmLog | cut -d : -f 1)

## Get snakemake jobids for jobs with Errors for splitFASTQ & alignFASTQ (jobs may appear multiple times if multiple errors occurred)
splitErrors=$(head -n $endSplit $slurmLog | grep 'Error' | grep -o 'jobid: [0-9]*' | sed 's/jobid: //')
alignErrors=$(tail -n +$endSplit $slurmLog | grep 'Error' | grep -o 'jobid: [0-9]*' | sed 's/jobid: //')

## Get snakemake jobids for Error jobs that finished
splitFinished=$(for i in $splitErrors; do grep -x "Finished job $i." <(head -n $endSplit $slurmLog) | sed 's/Finished job \([0-9]*\)\./\1/'; done)
alignFinished=$(for i in $alignErrors; do grep -x "Finished job $i." <(tail -n +$endSplit $slurmLog) | sed 's/Finished job \([0-9]*\)\./\1/'; done)

## Check that all Error jobs eventually finished
if [[ $splitErrors == $splitFinished ]]
then
    echo 'splitFASTQ jobs successful!'
else
    echo "Something isn't right with splitFASTQ jobs..."
fi

if [[ $alignErrors == $alignFinished ]]
then
    echo 'alignFASTQ jobs successful!'
else
    echo "Something isn't right with alignFASTQ jobs..."
fi