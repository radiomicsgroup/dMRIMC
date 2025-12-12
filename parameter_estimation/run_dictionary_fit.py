import os
import subprocess
from tqdm import tqdm
import sys

substrates = [
    "sub1",
    "sub2",
    "sub3",
    "sub4",
    "sub5",
    "sub6",
    "sub7",
    "sub8",
    "sub9",
    "sub10",
    "sub11",
    "sub12",
    "sub13",
    "sub14",
    "sub15",
    "sub16",
    "sub17",
    "sub18",
]

seqs = ["PGSEin"]
SNR = 50
configs = [14]
ncpu = sys.argv[1]
reg = 0.0000
successes = []
failures = []
total_runs = len(seqs) * len(configs) * len(substrates)
for seq in seqs:
    for config in configs:
        for sub in substrates:
            print(f"{total_runs} items left")
            output_path = f"dict_fit/config_{config}/{seq}_{sub}_out"
            data_path = f"LeaveOneOut"
            if not os.path.exists(output_path):
                os.makedirs(output_path)

            cmd = f"python mri2micro_dictml.py --noise {data_path}/noise_maps/noise_LOO_SNR{SNR}_{seq}_{sub}_out.nii --reg 2,{reg} --sldim 0 --ncpu {ncpu} {data_path}/niftis/signals/{seq}_sigs_noisy_SNR{SNR}_{sub}_out.nii {data_path}/npy_files/signals/{seq}_{sub}_out_train_SIGNALS.npy {data_path}/npy_files/parameters/{seq}_{sub}_out_train_PARAMS.npy {output_path}/SNR_{SNR}"
            print(cmd)

            print(f"Current job > Config = {config} | SNR = {SNR} | File = {seq}_{sub}_out")
            try:
                subprocess.run(cmd, shell=True, check=True,
                               stdout=subprocess.DEVNULL,
                               stderr=subprocess.DEVNULL)
                succ_msg = f"{seq}_{sub}_out terminated successfully"
                print(succ_msg)
                successes.append(succ_msg)
            except:
                fail_mess = f"{seq}_{sub}_out failed"
                print(fail_mess)
                failures.append(fail_mess)
            total_runs -= 1
print(successes)
print(failures)