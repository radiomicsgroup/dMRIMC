import numpy as np
import glob as gb
import pickle as pk
from tqdm import tqdm


def vol_w_intra(substrate, sequence, dest_folder, protocol):
    """
    Import the signals of all cells, volume weight them, add them and save the 
    resulting array.

    """
    # Areas
    lfile = open("metrics_and_cell_data/cells_areas.bin", 'rb')
    areas = pk.load(lfile)
    lfile.close()

    # Complete signals for each diffusivity D0
    # Signals
    sig_dir = f"{protocol}/{sequence}_{substrate}_cells"
    dval = np.linspace(0.8, 3, 10)
    dval = np.round(dval, decimals=4)

    for val in dval:
        print(f"Processing dvals of {val}")
        signal_files = gb.glob(f"{sig_dir}/*{val}.signals.bin")
        # Import data
        generated_signals = []
        for i in tqdm(signal_files):
            # print(i)
            s = i.split("/")
            name = s[2].split("_")[0]
            # will be the same for all probably but lets keep it
            D0 = s[2].split("_")[1][1:4]
            hfile = open(i, 'rb')
            signal = pk.load(hfile)
            hfile.close()
            generated_signals.append([signal[0][0], name, float(D0)])

        # Volume weighting
        sub_vols = areas[f'{substrate}']
        # Total volume
        tot_vol = sum(sub_vols.values())

        weighted_signals = []
        for i in tqdm(generated_signals):
            vol_frac = sub_vols[f'{i[1]}.ply'] / tot_vol
            w = i[0] * vol_frac
            # Removed because there should not be any NaNs
            # if np.any(np.isnan(w)):
            #     w = np.nan_to_num(w)
            nw = [w, i[1], i[2]]
            weighted_signals.append(nw)

        # Total intracellular signal
        # print(generated_signals)
        total_intra_signal = np.zeros_like(generated_signals[0][0])
        for i in weighted_signals:
            total_intra_signal += i[0]
        # Save array
        hfile = open(
            f"{dest_folder}/{sequence}_{substrate}_{nw[2]}_total_intra_signal.bin", 'wb')
        pk.dump(total_intra_signal, hfile, protocol=5)
        hfile.close()
