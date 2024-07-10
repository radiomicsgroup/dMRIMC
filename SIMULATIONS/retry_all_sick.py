import os
import sys
from parallel_execution_functions import run_bash_commands_in_parallel_with_timeout
import subprocess

sim_path = sys.argv[1]
try:
    os.chdir(sim_path)
except:
    ValueError("Simulation path unreachable or doesn't exist")

sub_type = sys.argv[2]
cores = int(sys.argv[3])


if sub_type == "intra":
    # For intracellular simulations the time window is 60 minutes --> 3600s
    time = 18000
elif sub_type == "extra":
    # For extracellular simulations the time window is 5 days --> 432000s
    time = 432000
else:
    ValueError("Substrate type must be <intra> or <extra>")
# Make copies of the error files because they are overwritten every time the sims run
# The function has been updated to make new files but never hurts to be safe
subprocess.run(['cp', "errors.txt", "errors_COPY.txt"])
subprocess.run(['cp', "timed_out.txt", "timed_out_COPY.txt"])
subprocess.run(['cp', "sickly_files.txt", "sickly_files_COPY.txt"])


def compare_error_files(original, new):
    """
    After we retry the files check if the contents of the new error files 
    are the same with the old ones
    """
    print(f"Comparing {original} with {new}")
    with open(f'{original}', 'r') as file1:
        file1_contents = file1.read().splitlines()

    with open(f'{new}', 'r') as file2:
        file2_contents = file2.read().splitlines()
    file1_contents.sort()
    file2_contents.sort()

    if set(file1_contents) == set(file2_contents):
        print("Error file contents are identical")
    else:
        print("Error file contents are different")
        diffs = set(file1_contents) - set(file2_contents)
        print(f"These files do not appear in the new file: \n{diffs}")
    return


# Sims that got stuck
print("Retrying those that got stuck or got sick for some other reason")
sick = open(f"sickly_files_COPY.txt", 'rt')
data3 = sick.read().splitlines()
sick.close()
# If we have no sick files run only timed out files
if data3 == []:
    print("No sick files, look again")
else:
    data4 = []  # sick files with same path as timed out for comparison purposes
    for i in data3:
        i = i.replace("random_walks", "conf_files")
        i = i.replace("_0.traj", ".conf")
        data4.append(i)
    print("Running some remedies...")
    commands = []
    for i in data4:
        proc = ["../MC-DC_Simulator_mod", f"{i}"]
        commands.append(proc)
    run_bash_commands_in_parallel_with_timeout(commands, n_parallel=cores, t=time, retry=True)
    compare_error_files("sickly_files.txt", "sickly_files_COPY.txt")
