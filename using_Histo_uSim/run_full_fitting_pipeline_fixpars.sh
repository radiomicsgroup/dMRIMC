#!/bin/bash

# In this example we process a scan containing multiple b-values acquired at one, single diffusion time
#
# The scan is stored in the https://github.com/radiomicsgroup/dMRIMC/tree/main/using_Histo_uSim/zenodo_mouse_data folder
# where you'll find:
# dwi_denoise_unring_sphmean_difftimefixed.nii --> scan with a fixed diffusion time (delta = 12 ms, Delta = 16.5 ms)
# dwi_denoise_unring_sphmean_difftimefixed.bval --> array of b-values (s/mm2) corresponding to volumes in dwi_denoise_unring_sphmean_difftimefixed.nii
# dwi_denoise_unring_sphmean_difftimefixed.gdur --> array of grad. duration (ms) corresponding to volumes in dwi_denoise_unring_sphmean_difftimefixed.nii
# dwi_denoise_unring_sphmean_difftimefixed.gsep --> array of grad. separation (ms) corresponding to volumes in dwi_denoise_unring_sphmean_difftimefixed.nii

#### Pipeline
echo "Combining signal and parameter arrays"

python combine_arrays.py

echo "Finding closest protocol, getting the subsets and processing .nii files"

python get_closest_scheme.py \
    --protocol-name MOUSE_BREAST_EXVIVO_FIXDIFFTIME \
    --ref-signal all_signals_all_substrates.npy \
    --dwi zenodo_mouse_data/dwi_denoise_unring_sphmean_difftimefixed.nii \
    --bval zenodo_mouse_data/dwi_denoise_unring_sphmean_difftimefixed.bval \
    --gdur zenodo_mouse_data/dwi_denoise_unring_sphmean_difftimefixed.gdur \
    --gsep zenodo_mouse_data/dwi_denoise_unring_sphmean_difftimefixed.gsep \
    --bval-threshold 20 \
    --vasc-threshold 250 \
    --noise zenodo_mouse_data/dwi_noise.nii \
    --show-plot

echo "Selecting parameter configuration: fin vCS_cyl D0in D0ex kappa"

python select_parameter_configuration.py --params fin vCS_cyl --output-folder MOUSE_BREAST_EXVIVO_FIXDIFFTIME

# Create dictionary with fixed kappa, for example kappa = 25 um/s and fixed D0in = 2.45
python fix_parameters.py --kappa 25 --D0in 2.45 --protocol-folder MOUSE_BREAST_EXVIVO_FIXDIFFTIME

echo "Fixed parameter array, the columns with the fixed values were removed"

echo "Making 'fitting' folder"

mkdir -v protocols/MOUSE_BREAST_EXVIVO_FIXDIFFTIME/fitting

echo "Running fitting script with fixed kappa = 25 μm/s and D0in = 2.45 μm2/ms"

python mri2micro_dictml.py \
    protocols/MOUSE_BREAST_EXVIVO_FIXDIFFTIME/dwi_normalized.nii \
    protocols/MOUSE_BREAST_EXVIVO_FIXDIFFTIME/kappa_D0in_fixed_sig.npy \
    protocols/MOUSE_BREAST_EXVIVO_FIXDIFFTIME/kappa_D0in_fixed_params.npy \
    --sldim 0 \
    --savg 3 \
    --ncpu 10 \
    --reg "2,0.0025" \
    --noise protocols/MOUSE_BREAST_EXVIVO_FIXDIFFTIME/dwi_noise_normalized.nii \
    --mask zenodo_mouse_data/dwi_mask_one_sample.nii \
    protocols/MOUSE_BREAST_EXVIVO_FIXDIFFTIME/fitting/Histo_uSim_difftimefixed
    
