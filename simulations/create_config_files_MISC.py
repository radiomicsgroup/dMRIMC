import numpy as np
import os
import glob as gb
import sys

"""
Create configuration files for the MCDC simulator

Example:
python create_config_files_MISC.py <sim_path> <substrate> <structure>

sim_path: folder that will hold the random walks
substrate: name of the substrate
structure: cells/misc/lumen etc
"""

# Path
sim_path = sys.argv[1]
try:
    os.chdir(sim_path)
except:
    ValueError("Simulation path unreachable or doesn't exist")
# Simulation folders
output_path = f"./conf_files"
randomwalks_path = f"./random_walks"
if not os.path.exists(output_path):
    os.makedirs(output_path)
if not os.path.exists(randomwalks_path):
    os.makedirs(randomwalks_path)

# Simulation parameters
Nvals = 1000
T = 3000
dur = 140 # ms
ply_scale = 1
dval = np.linspace(0.8, 3, 10)
dval = np.round(dval, decimals=4)
scheme_file = "../scheme_files/scheme.scheme"
substrate = sys.argv[2]
structure = sys.argv[3]
ply_directory = f"../../playgrounds/{substrate}/{structure}/plys"
ply_files = gb.glob(f"{ply_directory}/*.ply")

dim = []
with open(f"{ply_directory}/dimensions.txt") as f:
    for i in f:
        dim.append(i.split(','))
dims = np.array(dim)
sorted_dims = dims[dims[:, 0].argsort()]
ply_files = sorted(ply_files)
for j,i in enumerate(ply_files):
    xval = float(sorted_dims[j][1])
    yval = float(sorted_dims[j][2])
    cell_name = sorted_dims[j][0]
    for d in dval:
        with open(f'{output_path}/{sorted_dims[j][0]}_D{d}.conf', 'a') as file:
            l1 = f"N {Nvals}\n"
            l2 = f"T {T}\n"
            l3 = f"duration {dur}\n" # needs to be in seconds apparently
            l4 = f"diffusivity {d}e-6\n"
            l5 = f"exp_prefix {randomwalks_path}/{cell_name}_D{d}\n"
            l7 = f"scheme_file {scheme_file}\n"
            l8 = "\n"
            l9 = "write_txt 0\n"
            l10 = "write_bin 1\n"
            l11 = "write_traj_file 1\n"
            l12 = "\n"
            l13 = "<obstacle>\n"
            l14 = f"ply {i}\n"
            l15 = f"ply_scale {ply_scale}\n"
            l16 = "</obstacle>\n"
            # REMEMBER VOXELS ARE IN MILLIMETERS!!!!!!
            l17 = "\n"
            l18 = "<voxels>\n"
            # Adding a voxel big enough for any of the areas in x,y and 151μm in z
            # The paths have no negative values so start from 0
            l19 = "0.00 0.00 0.00\n"
            l20 = f"{xval + 0.1*xval} {yval + 0.1*yval} 0.151\n"
            l21 = "</voxels>\n"
            l22 = "\n"
            l23 = "<spawning_area>\n"
            # ALSO IN MILLIMETERS!!!!!!	    
            l24 = "0.000 0.000 0.06\n"
            # l25 = f"{xval + 0.01*xval} {yval + 0.01*yval} 0.09\n"
            l25 = f"{xval} {yval} 0.09\n"
            l26 = "</spawning_area>\n"
            l27 = "\n"
            l28 = "ini_walkers_pos intra\n"
            l29 = "\n"
            l30 = "num_process 1\n"
            l31 = "\n"
            l32 = "<END>\n"
            file.writelines([l1,l2,l3,l4,l5,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18,l19,l20,l21,l22,l23,l24,l25,l26,l27,l28,l29,l30,l31,l32])

# Simple counting tests for the times the IO is weird
print(f"{len(dval) * len(ply_files)} number of conf files expected")
print(f"{len(os.listdir(output_path))} is what exists in the target folder")
