import os
import subprocess
from tqdm import tqdm


substrates = ['substrate_1',
              'substrate_2',
              'substrate_3',
              'substrate_4',
              'substrate_5',
              'substrate_6',
              'substrate_7',
              'substrate_8',
              'substrate_9',
              'substrate_10',
              'substrate_11',
              'substrate_12',
              'substrate_13',
              'substrate_14',
              'substrate_15',
              'substrate_16',
              'substrate_17',
              'substrate_18',
              ]
seqs = ["PGSEin", "PGSEex", "TRSE"]
SNRs = ["50", "20"]
configurations = [9,13]
ncpu = 10
successes = []
failures = []
total_runs = len(seqs) * len(SNRs) * len(configurations) * len(substrates)
for seq in seqs:
    for SNR in SNRs:
        for config in configurations:
            for sub in substrates:
                print(f"{total_runs} items left")
                output_path = f"res/config_{config}/{seq}_{sub}_out"
                if not os.path.exists(output_path):
                    os.makedirs(output_path)
                cmd = f"python mri2micro_dictml.py --noise niftiis/noise_maps/noise_train_test_split_SNR{SNR}_{seq}_{sub}_out.nii --reg 2,0 --sldim 0 --ncpu {ncpu} niftiis/signals/{seq}_sigs_noisy_SNR_{SNR}_{sub}_out.nii signal_np_files/{seq}_{sub}_out_train.npy param_np_files/{seq}_param_arr_config{config}_{sub}_out_train.npy {output_path}/SNR_{SNR}"
                print(
                    f"Current job > Config = {config} | SNR = {SNR} | File = {seq}_{sub}_out")
                try:
                    subprocess.run(cmd, shell=True, check=True,
                                   stdout=subprocess.DEVNULL,
                                   stderr=subprocess.DEVNULL)
                    succ_mess = f"{seq}_{sub}_out terminated successfully"
                    print(succ_mess)
                    successes.append(succ_mess)

                except:
                    fail_mess = f"{seq}_{sub}_out failed"
                    print(fail_mess)
                    failures.append(fail_mess)
                total_runs -= 1


print(successes)
print(failures)
