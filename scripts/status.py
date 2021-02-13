#!/usr/bin/env python
import subprocess
import sys

jobid = sys.argv[1]

status = str(subprocess.check_output("sacct -j %s --format State --noheader | head -1 | awk '{print $1}'" % jobid, shell=True).strip())


## Adapted from: https://github.com/Snakemake-Profiles/slurm/blob/master/%7B%7Bcookiecutter.profile_name%7D%7D/slurm-status.py
if "BOOT_FAIL" in status:
  print("failed")
elif "OUT_OF_MEMORY" in status:
  print("failed")
elif "CANCELLED" in status:
  print("failed")
elif "COMPLETED" in status:
  print("success")
elif "DEADLINE" in status:
  print("failed")
elif "FAILED" in status:
  print("failed")
elif "NODE_FAIL" in status:
  print("failed")
elif "PREEMPTED" in status:
  print("failed")
elif "TIMEOUT" in status:
  print("failed")
# Unclear whether SUSPENDED should be treated as running or failed
elif "SUSPENDED" in status:
  print("failed")
else:
  print("running")