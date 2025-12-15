#!/bin/bash

echo "Combining signal and parameter arrays"

python combine_arrays.py

echo "Finding closest protocol, getting the subsets and processing .nii files"

python get_closest_scheme.py \
    --protocol-name MOUSE_BREAST_EXVIVO \
    --ref-signal all_signals_all_substrates.npy \
    --dwi zenodo_mouse_data/dwi_denoise_unring_sphmean.nii \
    --bval zenodo_mouse_data/dwi_denoise_unring_sphmean.bval \
    --gdur zenodo_mouse_data/dwi_denoise_unring_sphmean.gdur \
    --gsep zenodo_mouse_data/dwi_denoise_unring_sphmean.gsep \
    --bval-threshold 20 \
    --vasc-threshold 250 \
    --noise zenodo_mouse_data/dwi_noise.nii \
    --show-plot

echo "Selecting parameter configuration: fin vCS_cyl D0in D0ex kappa"

python select_parameter_configuration.py --params fin vCS_cyl

echo "Making 'fitting' folder"

mkdir -v fitting

echo "Running fitting script"

python mri2micro_dictml.py \
    dwi_normalized.nii \
    signal_arr_subset.npy \
    param_arr_subset.npy \
    --sldim 0 \
    --savg 3 \
    --ncpu 10 \
    --reg "2,0.0025" \
    --noise dwi_noise_normalized.nii \
    --mask zenodo_mouse_data/dwi_mask_one_sample.nii \
    fitting/Histo_uSim