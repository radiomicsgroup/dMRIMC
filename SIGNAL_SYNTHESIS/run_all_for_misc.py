import sys
import os
from MRI_functions import synthesis_misc_PGSE

"""
Run the signal synthesis for all substrates, both for intra and extra signals

Misc here stands for misc protocol, use the name of the 
sequence as it appears in the folder, for example:

run_all_for_misc.py <SEQUENCE_NAME> EXTRA
"""

substrates = os.listdir("../playgrounds/")
seq = sys.argv[1]
substrate_type = sys.argv[2] # <cells> or <EXTRA>
if substrate_type == "cells":
    Nspins = 1000
elif substrate_type == "EXTRA":
    Nspins = 10000
else:
    raise ValueError("substrate type must be <cells> or <EXTRA>")

for i in substrates:
    print(f"working on {i} for {seq}")
    synthesis_misc_PGSE(Nspins=Nspins, Ntime=3000, Tdur=0.140, substrate=i, substrate_type=substrate_type, snr=20, sequence=seq)
