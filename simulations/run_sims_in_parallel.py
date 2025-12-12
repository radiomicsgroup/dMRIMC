import os
from collections import deque
import sys
from parallel_execution_functions import run_bash_commands_in_parallel_with_timeout
import subprocess

"""
python run_sims_in_parallel.py <sub_name> <number of cores>
"""

# Path
sim_path = f"simulations/{sys.argv[1]}"
try:
    os.chdir(sim_path)
except:
    ValueError("Simulation path unreachable or doesn't exist")

files = os.listdir("conf_files")
if len(files) == 0:
    raise ValueError("There are no configuration files in this folder")

# Clean files from previous runs
subprocess.run([f"rm -rv ./random_walks/*"],shell=True)

# The commands are lists of the form [binary, target file]
commands = []
for i in files:
    proc = ["../../MC-DC_Simulator", f"conf_files/{i}"]
    commands.append(proc)

cores = int(sys.argv[2])
# Time window is 6 hours --> 21600s
run_bash_commands_in_parallel_with_timeout(commands, n_parallel=cores, t=21600)

# Print total errors and timeouts
errs = open(f"errors.txt")
tims = open(f"timed_out.txt")
derrs = errs.read().splitlines()
if derrs != []:
    print(f"{len(derrs)} failed")
else:
    print("No errors")
dtims = tims.read().splitlines()
if dtims != []:
    print(f"{len(dtims)} timed out")
else:
    print("No time outs")

