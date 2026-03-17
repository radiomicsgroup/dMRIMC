import numpy as np
import argparse

# Load combined parameter array
existing_array = np.load("all_parameters_all_substrates.npy")

# Load cellularity column from CommBio table 1
cellularity_column = np.load('reference_param_arrays/cellularity_column_from_CommBio.npy')

# Transform cellularity to log space
cellularity_column = np.log10(cellularity_column)

# Add the cellularity column
param_arr = np.hstack([existing_array, cellularity_column])


# Parameters are stored in this order:
# fin, mCS, varCS, skewCS, vCS_sph, vCS_cyl, Din, Dex, kappa, cellularity
param_index = {
    "fin":0,
    "mCS":1,
    "varCS":2,
    "skewCS":3,
    "vCS_sph":4,
    "vCS_cyl":5,
    "Din":6,
    "Dex":7,
    "kappa":8,
    "cellularity":9,
    }

# Argparse setup
parser = argparse.ArgumentParser(description="Select a parameter configuration for fitting")
parser.add_argument("--params", type=str, required=True, 
                    choices=["fin", "mCS", "varCS", "skewCS", "vCS_sph", 
                             "vCS_cyl", "Din", "Dex", "kappa", "cellularity"],
                    nargs="+",
                    help="Parameter names to keep")
parser.add_argument("--output-folder", type=str, help="Output folder to store the processed arrays")
args = parser.parse_args()

out_folder = args.output_folder

p_columns = args.params
# Add Din, Dex, kappa only if not already present
for p in ["Din", "Dex", "kappa"]:
    if p not in p_columns:
        p_columns.append(p)
cols_to_keep = [param_index[p] for p in p_columns]

param_arr_filtered = param_arr[:, cols_to_keep]

np.save(f"protocols/{out_folder}/param_arr_subset.npy", param_arr_filtered)
print(f"Saved parameter array subset at protocols/{out_folder}/param_arr_subset.npy")