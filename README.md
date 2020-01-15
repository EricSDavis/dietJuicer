# dietJuicer

Like `juicer`, but without all those calories

`dietJuicer` is a lighter-weight, HPC flexible version of juicer written with snakemake.


## Quickstart

The following instructions are tailored for the UNC longleaf cluster. For portability, change the indicated sections to match your cluster specifications.

1. Clone workflow into working directory
    ```
    git clone path/to/git/repo
    
    ```

## Setup & Dependencies

* python 3 package `psutil` for benchmarking
    ```
    module load python/3.6.6
    pip3.6 install --user psutil
    ```
