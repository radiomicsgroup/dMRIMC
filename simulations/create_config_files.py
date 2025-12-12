import numpy as np
import os
import glob as gb

# Simulation parameters
Nspins = 20000
T = 2370
dur = 110  # ms
ply_scale = 1
dval = np.linspace(0.8, 3, 5)
dval = np.round(dval, decimals=4)
permeability_value_range = np.linspace(0,40,9) # um/s
ply_files = gb.glob("ply_files/*.ply")

for sub in ply_files:
    # remove first part
    sub = sub.split("/")[1]
    substrate = sub[:-4] # get rid of ".ply"
    sim_path = f"simulations/{substrate}"
    if not os.path.exists(sim_path):
        os.makedirs(sim_path)
    # Simulation folders
    output_path = f"{sim_path}/conf_files"
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    randomwalks_path = f"{sim_path}/random_walks"
    if not os.path.exists(randomwalks_path):
        os.makedirs(randomwalks_path)

    # Substrate properties
    with open(f'ply_files/dims/{substrate}.txt', 'r') as file:
        properties = file.readlines()

    sub_length_x = float(properties[0].split(",")[1]) # mm
    sub_length_y = float(properties[0].split(",")[2]) # mm
    # Buffer between the substrate and the voxel in mm
    voxel_buffer_x = 0.1 * sub_length_x
    voxel_buffer_y = 0.1 * sub_length_y

    ply = substrate
    prefix = f"./random_walks/{substrate}"

    # print(ply)
    # print(prefix)
    # print(output_path)

    for kappa in permeability_value_range:
        for dvalin in dval:
            for dvalex in dval:
                # Distance an extracellular spin can traverse during the simulation
                spin_movement = np.sqrt(2 * dvalex * dur)
                spin_movement = spin_movement * 1e-3  # make it mm
                # Length of spawning area (half of the spawning area cube side)
                a = (sub_length_x / 2) - spin_movement

                current_filename = f'{output_path}/{substrate}_Din{dvalin}_Dex{dvalex}_k{kappa}um_per_s.conf'
                with open(current_filename, 'a') as file:
                    l1 = f"N {Nspins}\n"
                    l2 = f"T {T}\n"
                    l3 = f"duration {dur}\n" # ms
                    l_diffintra = f"diff_intra {dvalin}e-6\n" # mm^2/ms
                    l_diffextra = f"diff_extra {dvalex}e-6\n" # mm^2/ms
                    l5 = f"exp_prefix {prefix}_Din{dvalin}_Dex{dvalex}_k{kappa}um_per_s\n"
                    l7 = f"scheme_file ../../scheme_files/scheme.scheme\n"
                    l8 = "\n"
                    l9 = "write_txt 0\n"
                    l10 = "write_bin 1\n"
                    l11 = "write_traj_file 1\n"
                    l12 = "\n"
                    l13 = "<obstacle>\n"
                    l14 = f"ply ../../ply_files/{ply}.ply\n"
                    l15 = f"ply_scale {ply_scale}\n"
                    lperm = f"ply_permeability {kappa}e-6\n" # make it m/s
                    l16 = "</obstacle>\n"
                    l17 = "\n"
                    l18 = "<voxels>\n" # mm
                    # The paths have no negative values so start from 0
                    l19 = f"0.00 0.00 0.00\n"
                    l20 = f"{sub_length_x + voxel_buffer_x} {sub_length_y + voxel_buffer_y} 0.150\n"
                    l21 = "</voxels>\n"
                    l22 = "\n"
                    l23 = "<spawning_area>\n" # mm
                    l24 = f"{spin_movement} {spin_movement} 0.06\n"
                    l25 = f"{sub_length_x - spin_movement} {sub_length_y - spin_movement} 0.09\n"
                    l26 = "</spawning_area>\n"
                    l27 = "\n"
                    # l28 = "ini_walkers_pos extra\n"
                    l29 = "\n"
                    l30 = "num_process 1\n"
                    l31 = "\n"
                    l32 = "<END>\n"
                    file.writelines([l1, l2, l3, l_diffintra, l_diffextra, l5, l7, l8, l9, l10, l11, l12, l13, l14, l15, lperm, l16,
                                     l17, l18, l19, l20, l21, l22, l23, l24, l25, l26, l27, l29, l30, l31, l32])

    # Simple counting test for the times the IO is weird
    print(f"{len(dval)**2 * len(permeability_value_range)} number of conf files expected")
    print(f"{len(os .listdir(output_path))} is what exists in the target folder")
