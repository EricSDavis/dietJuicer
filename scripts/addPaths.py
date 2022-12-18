# -----------------------------------------------------------------------------
# Example usage:
# python3 ./scripts/addPaths.py\
#  -d "../MOMA_techreps/output/"\ # path to dietJuicerCore "output" directory
#  -s "samplesheet.txt"\ # path to the samplesheet
#  -c 'Project' 'Cell_Type' 'Genotype' 'Bio_Rep' 'Tech_Rep'\ # column names
#  -o "new_samplesheet.txt" # output (use same as -s to overwrite)
# -----------------------------------------------------------------------------

import argparse
import os
import pandas as pd
import glob

# Instantiate parser
parser = argparse.ArgumentParser(
    description="Add paths of resulting dietJuicerCore run to a samplesheet.")

# Add arguments
parser.add_argument('-d', '--directory',
                    type=str,
                    dest='dir',
                    action='store',
                    required=True,
                    help="Path to dietJuicerCore output directory.")

parser.add_argument('-s', '--samplesheet',
                    type=str,
                    dest='samplesheet',
                    action='store',
                    required=True,
                    help='Add paths to this samplesheet.')

parser.add_argument('-c', '--columnNames',
                    type=str,
                    nargs='+',
                    dest='columnNames',
                    action='store',
                    required=True,
                    help='A space-separated list of column names from the\
                            samplesheet that define the directory names used for\
                            dietJuicerCore output directories.\
                            This list can be copied from the "groupBy"\
                            parameter in config/config.yaml.\
                            Example: -c "Project" "Cell_Type" "Genotype" "Bio_Rep" "Tech_Rep"')

parser.add_argument('-o', '--output',
                    type=str,
                    dest='output',
                    action='store',
                    default='new_samplesheet.txt',
                    help='File name for output samplesheet.\
                            Default is "new_samplesheet.txt",\
                            use the same as "-s" to overwrite.')

# Parse args
args = parser.parse_args()
args.dir = os.path.abspath(args.dir)  # convert to absolute path

# Ensure file paths exist
if not(os.path.exists(args.dir)):
    parser.error("dietJuicerCore output directory does not exist.")

if not(os.path.exists(args.samplesheet)):
    parser.error("samplesheet does not exist.")

# Read samplesheet & convert columns to strings
samples = pd.read_csv(args.samplesheet, delimiter='\t')
samples = samples.astype(str)

# columnNames should be a subset of samplesheet column names
if not(set(args.columnNames).issubset(set(samples.columns))):
    parser.error("Some column names do not exist in the samplesheet.")

# Get output directory names
dirs = os.listdir(args.dir)
dirs.remove('logs_slurm')

# For each dietJuicerCore output directory...
for d in dirs:
    # Define the paths for merged_nodups, inter, and inter30 files
    mergedNoDups = glob.glob(os.path.join(
        args.dir, d, "*_dedup_merged_nodups.txt.gz"))
    inter = glob.glob(os.path.join(
        args.dir, d, "*_inter.txt"))
    inter30 = glob.glob(os.path.join(
        args.dir, d, "*_inter_30.txt"))

    # Create dictionary between column names and values
    filt = dict(zip(args.columnNames, d.split("_")))

    # Identify matching rows
    matchingRows = (samples[list(filt)] == pd.Series(filt)).all(axis=1)

    # Add paths to samplesheet
    samples.loc[matchingRows, 'merged_nodups'] = mergedNoDups
    samples.loc[matchingRows, 'inter'] = inter
    samples.loc[matchingRows, 'inter30'] = inter30

# Write results to a tsv file
samples.to_csv(args.output, sep='\t', index=False)
