import os
from collections import deque
import sys
from parallel_execution_functions import run_bash_commands_in_parallel_with_timeout


"""
python run_sims_in_parallel.py <folder> <substrate type> <number of cores>
substrate type = intra/extra
"""

# Path
sim_path = sys.argv[1]
try:
    os.chdir(sim_path)
except:
    ValueError("Simulation path unreachable or doesn't exist")

files = os.listdir("conf_files")
if len(files) == 0:
    raise ValueError("There are no configuration files in this folder")

# The commands are lists of the form [binary, target file]
commands = []
for i in files:
    proc = ["../MC-DC_Simulator_mod", f"conf_files/{i}"]
    commands.append(proc)

sub_type = sys.argv[2]
cores = int(sys.argv[3])
if sub_type == "intra":
    # For intracellular simulations the time window is 60 minutes --> 3600s
    run_bash_commands_in_parallel_with_timeout(commands, n_parallel=cores, t=3600)
elif sub_type == "extra":
    # For extracellular simulations the time window is 5 days --> 432000s
    run_bash_commands_in_parallel_with_timeout(commands, n_parallel=cores, t=432000)
else:
    ValueError("Substrate type must be <intra> or <extra>")
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

