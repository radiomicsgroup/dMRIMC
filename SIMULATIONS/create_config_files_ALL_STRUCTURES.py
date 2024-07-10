import numpy as np
import os
import sys

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
dval = np.linspace(0.8, 3, 5)
dval = np.round(dval, decimals=4)
scheme_file = "../scheme_files/scheme.scheme"
substrate = sys.argv[2]
with open(f'../../playgrounds/{substrate}/dimensions.txt', 'r') as dimfile:
    dimfile_contents = dimfile.read().splitlines()

sub_length_x = float(dimfile_contents[0].split(",")[1]) # mm
sub_length_y = float(dimfile_contents[0].split(",")[2])

# Buffer between the substrate and the voxel in mm
voxel_buffer_x = 0.1 * sub_length_x
voxel_buffer_y = 0.1 * sub_length_y
# Spawning area adjustment
# The spawning area will be configured to be determined by the value
# of d. It should be contained in the substrate, not the other way
# around. So the spawning area will be dim(x/y) of substrate - the distance
# a spin can travel in the simulation duration.
ply = f"../../playgrounds/{substrate}/{substrate}_ALL_STRUCTURES.ply"
prefix = f"{randomwalks_path}/{substrate}_ALL_STRUCTURES"

for d in dval:
    # Distance a spin can traverse during the simulation 
    spin_movement = np.sqrt(2 * d * dur)
    spin_movement = spin_movement * 1e-3 # make it mm
    with open(f'{output_path}/{substrate}_ALL_STRUCTURES_D{d}.conf', 'a') as file:
        l1 = f"N {Nvals}\n"
        l2 = f"T {T}\n"
        l3 = f"duration {dur}\n"
        l4 = f"diffusivity {d}e-6\n"
        l5 = f"exp_prefix {randomwalks_path}/{substrate}_ALL_STRUCTURES_D{d}\n"
        l7 = f"scheme_file {scheme_file}\n"
        l8 = "\n"
        l9 = "write_txt 0\n"
        l10 = "write_bin 1\n"
        l11 = "write_traj_file 1\n"
        l12 = "\n"
        l13 = "<obstacle>\n"
        l14 = f"ply {ply}\n"
        l15 = f"ply_scale {ply_scale}\n"
        l16 = "</obstacle>\n"
        # REMEMBER VOXELS ARE IN MILLIMETERS!!!!!!
        l17 = "\n"
        l18 = "<voxels>\n"
        # Adding a voxel big enough for any of the areas in x,y and 151μm in z
        # The paths have no negative values so start from 0
        l19 = f"0.00 0.00 0.00\n"
        l20 = f"{sub_length_x + voxel_buffer_x} {sub_length_y + voxel_buffer_y} 0.150\n"
        l21 = "</voxels>\n"
        l22 = "\n"
        l23 = "<spawning_area>\n"
        # ALSO IN MILLIMETERS!!!!!!
        l24 = f"{spin_movement} {spin_movement} 0.06\n"
        l25 = f"{sub_length_x - spin_movement} {sub_length_y - spin_movement} 0.09\n"
        l26 = "</spawning_area>\n"
        l27 = "\n"
        l28 = "ini_walkers_pos extra\n"
        l29 = "\n"
        l30 = "num_process 1\n"
        l31 = "\n"
        l32 = "<END>\n"
        file.writelines([l1,l2,l3,l4,l5,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18,l19,l20,l21,l22,l23,l24,l25,l26,l27,l28,l29,l30,l31,l32])

# Simple counting test for the times the IO is weird
print(f"{len(dval)} number of conf files expected")
print(f"{len(os.listdir(output_path))} is what exists in the target folder")
