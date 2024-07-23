import os
import glob as gb
import re
import sys

def sim_tests(sim_type):
    files = os.listdir("conf_files")
    if len(files) == 0:
        raise ValueError("There are no configuration files in this folder")
    # Tests
    walks_per_file = len(gb.glob("random_walks/*.traj"))
    total_files = len(os.listdir('random_walks'))

    # First of all check for files that are not the correct number of bytes
    # With the current (11/5/2023) sim parameters the size for 1000 spins is
    # 36012000 but you need to input the number that you are expecting
    print("\n#######################")
    print("Checking for sick files")
    print("#######################")
    sfiles = open("sickly_files.txt", 'wt')
    if sim_type == "intra":
        byte_size = 36012000
    elif sim_type == "extra":    
        byte_size = 360120000
    else:
        raise ValueError("Simulation type must be <intra> or <extra>")
    traj_files = gb.glob("random_walks/*.traj")
    cnt_sick = 0
    for i in traj_files:
        file_size = os.path.getsize(i)
        if file_size < (byte_size):
            cnt_sick += 1
            sfiles.write(f"{i}\n")
            os.rename(i,f"{i}"+"_SICK")
    sfiles.close()
    print(f"There were {cnt_sick} sick files here")

    if walks_per_file == len(files):
        print("All good mate")
    elif walks_per_file != len(files):
        print(f"\nSomething's fishy, the number of walks is {walks_per_file} \nand we were expecting {len(files)}")
        print(f"Total number of files: {total_files}/{len(files)*5}")
        
        print("\nChecking to see if the trajectory files are messed up...")
        # Removing all files that are repeats, i.e have "rep" in their name
        # regular expression matching
        print("################################")
        print("Looking for trajectory rep files")
        print("################################")
        traj_files = gb.glob("random_walks/*.traj")
        cnt_reps = 0
        pattern = r"\b_rep_\b"
        for i in traj_files:
            match = re.search(pattern, i)
            # This part should be finding nothing, as all these should have been deleted
            # as sick files above
            if match and (file_size < (byte_size)):
                os.rename(i,f"{i}"+"_REP")
                cnt_reps += 1
            # This part will find those that are not sick but are still repeats
            elif match and not (file_size < (byte_size)):
                print("Rep file found, full size though, check manually")
                cnt_reps += 1
        print(f"Found {cnt_reps} rep files of full size\n")

        # TRAJECTORY FILES
        # Re-read the directory because changes might have been made from above
        traj_files = gb.glob("random_walks/*.traj")
        extras = open("extra_traj_files.txt", 'wt')
        missing = open("missing_traj_files.txt", 'wt')
        # Remove the surrounding text
        tfiles = sorted([i[13:] for i in traj_files])
        ffiles = sorted([i[:-5] for i in files])
        string = "_0.traj"
        # We expect to find only _0.traj files if all was done ok so we add
        # that string at the end of the filenames to cross reference them
        ffiles = [f"{element}{string}" for element in ffiles]

        pasta = set(tfiles) - set(ffiles)
        # print(pasta)
        if pasta == set():
            print("All trajectory files are there")
            print("The discrepancy stems from the other files it seems so we don't care")
        elif set(ffiles) < set(tfiles):
            print("There are excess files, probably some sims have run more than once")
            extra_files = set(tfiles) - set(ffiles)
            extras.write("".join([f"{i}\n" for i in extra_files]))
            # Delete the extra files
            print("\n############################")
            print("Deleting any redundant files")
            print("############################")
            # print(extra_files)
            cnt_redundant = 0
            for i in traj_files:
                if i[13:] in extra_files:
                    os.rename(i,f"{i}"+"_EXTRA")
                    cnt_redundant += 1
            print(f"Found {cnt_redundant} redundant files")
            # Check that it actually worked
            traj_files2 = gb.glob("random_walks/*.traj")
            tfiles2 = sorted([i[13:] for i in traj_files2])
            print(tfiles2)
            print(ffiles)
            pasta2 = set(tfiles2) == set(ffiles)
            print(f"Are things better now? --> {pasta2==set()}")
        elif len(tfiles) < len(ffiles):
            print("Trajectory files are missing")
            missing_files_trajectory = set(ffiles) - set(tfiles)
            missing.write("conf_files/".join([f"{i}\n" for i in missing_files_trajectory]))
        extras.close()
        missing.close()