import numpy as np
import argparse

D0in_D0ex_values = np.round(np.linspace(0.8, 3, 5), decimals=4)
kappa_values = np.linspace(0, 40, 9)

parser = argparse.ArgumentParser(description="Fix kappa, D0in, D0ex or a combination of them. ")
parser.add_argument("--kappa", type=float, choices=kappa_values,
                    help=f"Value for kappa (allowed: {kappa_values})")
parser.add_argument("--D0in", type=float, choices=D0in_D0ex_values,
                    help=f"Value for D0in (allowed: {D0in_D0ex_values})")
parser.add_argument("--D0ex", type=int, choices=D0in_D0ex_values,
                    help=f"Value for D0ex (allowed: {D0in_D0ex_values})")
parser.add_argument("--protocol-folder", type=str, help="Folder containing the signal and parameter array subsets for this protocol, also serves as ouput folder to store the processed arrays")
args = parser.parse_args()

if not any([args.kappa, args.D0in, args.D0ex]):
    parser.error("At least one of --kappa, --D0in, or --D0ex must be specified")

protocol_folder = args.protocol_folder

kappa_col = -1
D0ex_col = -2
D0in_col = -3

sig_arr = np.load(f"protocols/{protocol_folder}/signal_arr_subset.npy")
param_arr = np.load(f"protocols/{protocol_folder}/param_arr_subset.npy")

# dictionary of provided arguments and their column indices
provided_args = {
    "kappa": (args.kappa, kappa_col),
    "D0in": (args.D0in, D0in_col),
    "D0ex": (args.D0ex, D0ex_col)
}

# filter out None values (arguments not provided)
provided_args = {k: v for k, v in provided_args.items() if v[0] is not None}

# print what is being fixed
print("Fixing:")
for arg, (val, col) in provided_args.items():
    print(f"{arg}: {val}")

# initialize mask
mask = np.ones(len(param_arr), dtype=bool)

# update mask for each provided argument
for val, col in provided_args.values():
    mask &= (param_arr[:, col] == val)

# apply the combined mask to both arrays
param_arr_fixed = param_arr[mask]
sig_arr_fixed = sig_arr[mask]

# remove columns corresponding to fixed parameters from the
# parameter array
cols_to_remove = [col for _, col in provided_args.values()]
param_arr_cleaned = np.delete(param_arr_fixed, cols_to_remove, axis=1)

# make the name for the output files
name = "_".join(provided_args.keys())
output_param_path = f"protocols/{protocol_folder}/{name}_fixed_params.npy"
output_sig_path = f"protocols/{protocol_folder}/{name}_fixed_sig.npy"

# Save the new arrays
np.save(output_param_path, param_arr_cleaned)
np.save(output_sig_path, sig_arr_fixed)

print()
print(f"Saved fixed parameter array to: {output_param_path}")
print(f"Saved fixed signal array to: {output_sig_path}")