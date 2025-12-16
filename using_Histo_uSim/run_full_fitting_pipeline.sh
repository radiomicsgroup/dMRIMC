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

python select_parameter_configuration.py --params fin vCS_cyl --output-folder MOUSE_BREAST_EXVIVO

echo "Making 'fitting' folder"

mkdir -v protocols/MOUSE_BREAST_EXVIVO/fitting

echo "Running fitting script"

python mri2micro_dictml.py \
    protocols/MOUSE_BREAST_EXVIVO/dwi_normalized.nii \
    protocols/MOUSE_BREAST_EXVIVO/signal_arr_subset.npy \
    protocols/MOUSE_BREAST_EXVIVO/param_arr_subset.npy \
    --sldim 0 \
    --savg 3 \
    --ncpu 10 \
    --reg "2,0.0025" \
    --noise protocols/MOUSE_BREAST_EXVIVO/dwi_noise_normalized.nii \
    --mask zenodo_mouse_data/dwi_mask_one_sample.nii \
    protocols/MOUSE_BREAST_EXVIVO/fitting/Histo_uSim