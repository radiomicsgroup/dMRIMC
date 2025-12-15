#!/bin/bash

# First we ran things like before to create the parameter and signal arrays
# and to process the dwi file
echo "Combining signal and parameter arrays"

python combine_arrays.py

echo "Finding closest protocol, getting the subsets and processing .nii files"

python get_closest_scheme.py \
    --protocol-name MOUSE_BREAST_EXVIVO \
    --ref-signal all_signals_all_substrates.npy \
    --dwi zenodo_mouse_data/scans/breast/dwi_denoise_unring_sphmean.nii \
    --bval zenodo_mouse_data/scans/breast/dwi_denoise_unring_sphmean.bval \
    --gdur zenodo_mouse_data/scans/breast/dwi_denoise_unring_sphmean.gdur \
    --gsep zenodo_mouse_data/scans/breast/dwi_denoise_unring_sphmean.gsep \
    --bval-threshold 20 \
    --vasc-threshold 250 \
    --noise zenodo_mouse_data/scans/breast/dwi_noise.nii \
    --show-plot

echo "Selecting parameter configuration: fin vCS_cyl D0in D0ex kappa"

python select_parameter_configuration.py --params fin vCS_cyl

# As an example, we create a small sub-protocol of diffusion MRI measurements at fixed diffusion time
echo "Selecting the first diffusion time from protocols/MOUSE_BREAST_EXVIVO/scheme as an example"
echo "Using signal_arr_subset_one_diff_time.npy"
echo "Grabbing the appropriate volumes from the normalized scan"

fslroi dwi_normalized.nii dwi_normalized_first_diff_time.nii 0 3

echo "0.00 500.00 2100.00" > protocols/MOUSE_BREAST_EXVIVO_ONE_DIFF_TIME/scheme_fixed_diffusion time.scheme
echo "11.00 11.00 11.00" > protocols/MOUSE_BREAST_EXVIVO_ONE_DIFF_TIME/scheme_fixed_diffusion time.scheme
echo "16.00 16.00 16.00" > protocols/MOUSE_BREAST_EXVIVO_ONE_DIFF_TIME/scheme_fixed_diffusion time.scheme

# Create dictionary with fixed kappa, for example kappa = 25 um/s and fixed D0in = 2.45
python fix_parameters.py --kappa 25 --D0in 2.45

echo "Fixed parameter array, the columns with the fixed values were removed"

echo "Making 'fitting' folder"

mkdir -v fitting

echo "Running fitting script with fixed kappa = 25 μm/s and D0in = 2.45 μm2/ms"

python mri2micro_dictml.py \
    dwi_normalized_first_diff_time.nii \
    kappa_D0in_fixed_sig.npy \
    kappa_D0in_fixed_params.npy \
    --sldim 0 \
    --savg 3 \
    --ncpu 10 \
    --reg "2,0.0025" \
    --noise dwi_noise_normalized.nii \
    --mask zenodo_mouse_data/scans/breast/dwi_mask_one_sample.nii \
    fitting/Histo_uSim