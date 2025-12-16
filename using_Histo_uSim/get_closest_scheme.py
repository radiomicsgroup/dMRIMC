import numpy as np
import os
import nibabel as nib
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter
import argparse

def fancier_print(msg):
    print(25 * "-")
    print(msg)
    print(25 * "-")
    return

parser = argparse.ArgumentParser(description="Get closest signal entries to target protocol, process and normalize DWI scans")
parser.add_argument("--protocol-name", required=True, help="Output protocol name")
parser.add_argument("--dwi", required=True, help="Input DWI NIfTI file")
parser.add_argument("--bval", required=True, help="Target bval file")
parser.add_argument("--gdur", required=True, help="Target gdur file")
parser.add_argument("--gsep", required=True, help="Target gsep file")
parser.add_argument("--ref-signal", required=True, default="all_signals_all_substrates.npy", help="Reference signal array")
parser.add_argument("--vasc-threshold", required=False, default="250", help="Remove the bvalues below this number to get rid of vascular signals")
parser.add_argument("--bval-threshold", required=False, default="20", help="bvalue threshold under which all volumes are considered b = 0")
parser.add_argument("--noise-map", required=False, help="OPTIONAL - NIfTI with the standard deviation of the noise. If supplied it will be normalized like the DWI file ")
parser.add_argument("--show-plot", action="store_true", help="OPTIONAL - Plot first row of final signal array, separated via diffusion time and/or gradient directions")
# action="store_true" makes including --show-plot --> args.show_plot = True, omitting it sets False
args = parser.parse_args()

# Load reference protocol
signal_arr = np.load(args.ref_signal)
ref_bvals = np.loadtxt(f"protocols/reference_protocol/ref.bval")
ref_delta = np.loadtxt(f"protocols/reference_protocol/ref.gdur")
ref_DELTA = np.loadtxt(f"protocols/reference_protocol/ref.gsep")

# Load target protocol
target_bvals = np.loadtxt(args.bval)
target_delta = np.loadtxt(args.gdur)
target_DELTA = np.loadtxt(args.gsep)

INPUT_PROTOCOL_NAME = args.protocol_name
os.makedirs(os.path.dirname(f"protocols/{INPUT_PROTOCOL_NAME}"), exist_ok=True)
#%%
dwi = nib.load(args.dwi)
data = dwi.get_fdata()
affine = dwi.affine
header = dwi.header
#%%
max_ref_bval = np.max(ref_bvals)
if np.any(target_bvals > max_ref_bval):
    print("Warning!")
    print(f"Max bvalue in the target scheme is {np.max(target_bvals)} s/mm2")
    print(f"Max bvalue synthesized signals are at {max_ref_bval} s/mm2 ")
    print(f"Volumes with bvalues higher than that are removed")
    too_high_locs = np.where(target_bvals > max_ref_bval)[0]
    data = np.delete(data, too_high_locs, axis=3)
    target_bvals = np.delete(target_bvals, too_high_locs)
    target_delta = np.delete(target_delta, too_high_locs)
    target_DELTA = np.delete(target_DELTA, too_high_locs)

vasc_bval_limit = int(args.vasc_threshold)
if np.any((target_bvals < vasc_bval_limit) & (target_bvals > int(args.bval_threshold))):
    print()
    print("Warning!")
    print(f"bvalue under {vasc_bval_limit} s/mm2 in the target scheme")
    print(f"Min bvalue is 300 s/mm2 so that all effects of vascular flow are removed")
    print(f"Volumes with bvalues lower than that (expect for b = 0) are removed")
    too_low_locs = np.where((target_bvals <= vasc_bval_limit) & (target_bvals > 20))[0]
    data = np.delete(data, too_low_locs, axis=3)
    target_bvals = np.delete(target_bvals, too_low_locs)
    target_delta = np.delete(target_delta, too_low_locs)
    target_DELTA = np.delete(target_DELTA, too_low_locs)
#%%
# Level 1: Filter by closest b
if np.any(target_bvals > np.max(ref_bvals)):
    print("Warning!")
    print(f"Max bvalue in the target scheme is {np.max(target_bvals)}")
    print(f"Max bvalue for the synthesized signals is {np.max(ref_bvals)}")
b_min_inds = []
for b in target_bvals:
    diff1 = np.abs(ref_bvals - b).argmin()
    b_min_inds.append(diff1)
res_bvals = ref_bvals[b_min_inds]
# if there are more than one b=0 in succession, keep only the first and save
# the indices of the ones we removed so we can remove those volumes from the
# scans as well
zero_indices = np.where(res_bvals == 0)[0]
consecutive_zeros = np.where(np.diff(zero_indices) == 1)[0] + 1
if consecutive_zeros.size > 0:
    indices_to_remove = zero_indices[consecutive_zeros]
    print("Some consecutive b = 0 volumes were created due to bvalues under 300")
    print(f"Removed extra b = 0 indices at {indices_to_remove}") 
    
    # update bmin indices
    b_min_inds = [i for j, i in enumerate(b_min_inds) if j not in indices_to_remove]
    # update res_bvals
    res_bvals = np.delete(res_bvals, indices_to_remove)
#%%
# Level 2
delta1_min_inds = []
for d1 in target_delta:
    diff2 = np.abs(ref_delta - d1).argmin()
    delta1_min_inds.append(diff2)
res_delta = ref_delta[delta1_min_inds]

# remove the volumes we removed from the bvals
if consecutive_zeros.size > 0:
    res_delta = np.delete(res_delta, indices_to_remove)
    delta1_min_inds = [i for j, i in enumerate(delta1_min_inds) if j not in indices_to_remove]
#%%
# Level 3
DELTA_min_inds = []
for D in target_DELTA:
    diff4 = np.abs(ref_DELTA - D).argmin()
    DELTA_min_inds.append(diff4)
res_DELTA = ref_DELTA[DELTA_min_inds]

# remove the volumes we removed from the bvals
if consecutive_zeros.size > 0:
    res_DELTA = np.delete(res_DELTA, indices_to_remove)
    DELTA_min_inds = [i for j, i in enumerate(DELTA_min_inds) if j not in indices_to_remove]
#%%
print()
fancier_print("Resulting scheme")
print(res_bvals)
print(res_delta)
print(res_DELTA)
print()
fancier_print("Target scheme")
print(target_bvals)
print(target_delta)
print(target_DELTA)
#%%
if not os.path.exists(f"protocols/{INPUT_PROTOCOL_NAME}"):
    os.makedirs(f"protocols/{INPUT_PROTOCOL_NAME}")
# Save resulting scheme
# Print if we removed any consecutive b=0
res_bvals = res_bvals.reshape(1,-1)
res_delta = res_delta.reshape(1,-1)
res_DELTA = res_DELTA.reshape(1,-1)
np.savetxt(f"protocols/{INPUT_PROTOCOL_NAME}/closest.bval", res_bvals, fmt="%.2f", delimiter=" ")
np.savetxt(f"protocols/{INPUT_PROTOCOL_NAME}/closest.gdur", res_delta, fmt="%.2f", delimiter=" ")
np.savetxt(f"protocols/{INPUT_PROTOCOL_NAME}/closest.gsep", res_DELTA, fmt="%.2f", delimiter=" ")

scheme = np.concatenate((res_bvals, res_delta, res_DELTA))
np.savetxt(f"protocols/{INPUT_PROTOCOL_NAME}/closest.scheme", scheme, fmt="%.2f", delimiter=" ")

print()
print(f"Saved closest protocol at: protocols/{INPUT_PROTOCOL_NAME}")
#%%
# Nifti processing
# normalize by b = 0, assuming it is the first volume. If there are many b = 0
# volumes, normalize by mean b = 0
# also do the noise
# Noise
noise = nib.load(args.noise_map)
noise_data = noise.get_fdata()
noise_affine = noise.affine
noise_header = noise.header
print()
print(f"Normalizing the scan and noise, assuming that any volumes with b < {args.bval_threshold} s/mm2 are b = 0")
b0_vols = np.where(target_bvals <= int(args.bval_threshold))[0]
if b0_vols.size > 1:
    mean_b0 = np.mean(data[..., b0_vols], axis=-1)
    new_data = data/mean_b0[..., np.newaxis]
    new_noise_data = noise_data/mean_b0
else:
    new_data = data/data[..., b0_vols]
    new_noise_data = noise_data/b0_vols

# Save result
output_img = nib.Nifti1Image(new_data, affine, header)
nib.save(output_img, f"protocols/{INPUT_PROTOCOL_NAME}/dwi_normalized.nii")

output_noise = nib.Nifti1Image(new_noise_data, noise_affine, noise_header)
nib.save(output_noise, f"protocols/{INPUT_PROTOCOL_NAME}/dwi_noise_normalized.nii")

if consecutive_zeros.size > 0:
    new_data = np.delete(new_data, indices_to_remove, axis=3)
    print("Removing the consecutive b = 0 volumes from the scan as well")

    # Save result
    edited_img = nib.Nifti1Image(new_data, affine, header)
    nib.save(edited_img, f"protocols/{INPUT_PROTOCOL_NAME}/dwi_normalized.nii")
#%% Grab signal subset
# Find indices that match ALL four parameters
b_min_inds_arr = np.array(b_min_inds)
delta1_min_inds_arr = np.array(delta1_min_inds)
DELTA_min_inds_arr = np.array(DELTA_min_inds)

# Find the reference indices that match our final closest scheme
final_inds = []
for i in range(res_bvals.size):
    # Find where ALL parameters match in the reference arrays
    if res_bvals[:, i] != 0:
        candidates = np.where(
            (ref_bvals == res_bvals[:, i]) &
            (ref_delta == res_delta[:, i]) &
            (ref_DELTA == res_DELTA[:, i])
        )[0]
        # print(candidates)
        if len(candidates) > 0:
            final_inds.append(candidates[0])
        else:
            raise ValueError(f"No matching entry found for volume {i}")
    else:
        # when we have a b = 0 volume S = 1
        final_inds.append(0)

# Extract signal subset
signal_subset = signal_arr[:, final_inds]
np.save(f"protocols/{INPUT_PROTOCOL_NAME}/signal_arr_subset.npy", signal_subset)
#%%
def split_by_zeros(arr1, arr2):
    """Split arr2 into rows based on positions of zeros in arr1."""
    zero_indices = np.where(arr1 == 0)[0]

    # Add endpoint to capture last segment
    split_points = np.append(zero_indices, len(arr1))

    rows = []
    for i in range(len(split_points) - 1):
        start = split_points[i]
        end = split_points[i + 1]
        rows.append(arr2[start:end])

    return np.array(rows)

if args.show_plot:
    sigs = split_by_zeros(res_bvals[0,:], signal_subset[0])
    bs = split_by_zeros(res_bvals[0,:], res_bvals[0,:])
    unsorted_found = False
    for i, s in enumerate(sigs):
        if np.any(bs[i] != np.sort(bs[i])):
            if not unsorted_found:
                plt.xticks(range(len(bs[i])), bs[i].astype(int))
                plt.xlabel("bvalue (s/mm$^2$)")
                unsorted_found = True
            plt.plot(s, label=f"Diffusion time {i+1}")
        else:
            plt.plot(bs[i], s, label=f"Diffusion time {i+1}")
            ticks = bs[i]
            ax = plt.gca()
            ax.xaxis.set_major_locator(FixedLocator(ticks))
            plt.xlabel("b-values (s/mm$^2$)")

    plt.ylabel("Signal intensity")
    plt.title("Signal at 0 index")
    plt.legend()
    plt.show()