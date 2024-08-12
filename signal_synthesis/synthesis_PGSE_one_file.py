import numpy as np
import glob as gb
import pickle as pk
import time
import os
from tqdm import tqdm
from MRI_functions import construct_signal_PGSE

Ntime = 3000
Tdur = 0.140
snr = 20

substrate = "mouse_example"
substrate_type = "cells"

if substrate_type == "cells":
    Nspins = 1000
elif substrate_type == "EXTRA":
    Nspins = 10000
else:
    raise ValueError("substrate type must be <cells> or <EXTRA>")


# Simulation parameters
Nsteps = Ntime + 1  # number of simulation interations
# Random walks directory
if substrate_type == "cells":
    trajdir = f"../simulations/{substrate}_CELLS"
elif substrate_type == "EXTRA":
    trajdir = f"../simulations/{substrate}_EXTRA"
# Get all trajectory files
traj_files = sorted(gb.glob(f"{trajdir}/random_walks/*.traj"))
# Sequence parameters
sequence = "CUSTOM_PGSE"
outstring = f"CUSTOM_PGSE/{sequence}_{substrate}_{substrate_type}"
if not os.path.exists(outstring):
    os.makedirs(outstring)

sequence_params_path = f"./sequence_parameters/{sequence}"
# List of b-value in sec/mm2
bvalseq = np.loadtxt(f"{sequence_params_path}/custom.bval", delimiter=' ')
# List of delta1 in msec
delta1seq = np.loadtxt(f"{sequence_params_path}/custom.gdur1", delimiter=' ')
# List of delta2 in msec
delta2seq = np.loadtxt(f"{sequence_params_path}/custom.gdur2", delimiter=' ')
# List of Delta1,2 in msec
Delta12seq = np.loadtxt(f"{sequence_params_path}/custom.gsep12", delimiter=' ')

########## START OF PROCESSING ##########
print(f"Starting signal synthesis")
print(f"Signals to construct: {len(traj_files)}")
start = time.time()
# Log file for errors
errors = 0
f = open(f"{outstring}/errors.txt", 'at')
for j in (pbar := tqdm(traj_files, total=len(traj_files))):
    traj_file = j.replace(trajdir, '')
    pbar.set_description(f"Processing trajectories for file {traj_file}")
    p = traj_file.replace(f"{trajdir}/", '')
    p = p.replace("_0.traj", '')
    p = p.replace("/random_walks", '')
    # Parameter list
    processinglist = [snr, bvalseq, delta1seq,
                      delta2seq, Delta12seq, Nspins, Nsteps, Tdur, j]
    # Compute MRI signals
    fitresults = construct_signal_PGSE(processinglist)
    outlist = [fitresults, p]
    hfile = open(f'{outstring}/{p}.signals.bin', 'wb')
    pk.dump(outlist, hfile, protocol=5)
    hfile.close()
f.close()
end = time.time()
print(f"Done")
print(f"Total time: {(end-start):.2f}s")
print(f"Total errors: {errors}")