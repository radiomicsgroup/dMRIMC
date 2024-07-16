import numpy as np
import glob as gb
import pickle as pk
from tqdm import tqdm


def vol_w_extra(substrate, sequence, dest_folder, protocol):
    """
    Import the signals of all cells, volume weight them, add them and save the 
    resulting array.
    """
    # Complete signals for each diffusivity D0
    # Signals
    sig_dir = f"{protocol}/{sequence}_{substrate}_EXTRA"
    dval = np.linspace(0.8, 3, 10)
    dval = np.round(dval, decimals=4)
    for val in dval:
        print(f"Processing dvals of {val}")
        signal_files_EXTRA = gb.glob(f"{sig_dir}/*{val}.signals.bin")

        generated_signals_EXTRA = []
        for i in tqdm(signal_files_EXTRA):
            s = i.split("/")
            name = substrate + "_EXTRA"
            # will be the same for all probably but lets keep it
            D0 = s[2].split("_")[-1][1:4]
            hfile = open(i, 'rb')
            signal = pk.load(hfile)
            hfile.close()
            generated_signals_EXTRA.append([signal[0][0], name, float(D0)])

        # Volume weighting
        weighted_signals = []
        weighted_signals = generated_signals_EXTRA[0]

        # Total extracellular signal
        total_extra_signal = np.zeros_like(generated_signals_EXTRA[0][0])
        total_extra_signal = weighted_signals[0]
        # Save array
        hfile = open(
            f"{dest_folder}/{sequence}_{substrate}_{weighted_signals[2]}_total_extra_signal.bin", 'wb')
        pk.dump(total_extra_signal, hfile, protocol=5)
        hfile.close()
