
## Histology-informed microstructural diffusion simulations for MRI cancer characterisation — the Histo-μSim framework

<div align="center">
  <img src="https://github.com/radiomicsgroup/dMRIMC/blob/main/imgs/commbio25.png" alt="commbio" width="auto" height="auto">
</div>

If you find dMRIMC useful, please cite our article:

Athanasios Grigoriou, et al. **"Histology-informed microstructural diffusion simulations for MRI cancer characterisation — the Histo-μSim framework"**. Communications Biology 2025, 8: 1695, doi: [10.1101/2024.07.15.24310280](https://doi.org/10.1038/s42003-025-09096-3).

This repository was developed by Athanasios Grigoriou (<agrigoriou@vhio.net>) and Francesco Grussu (<fgrussu@vhio.net>). **The project that gave rise to these results received the support of a fellowship from ”la Caixa” Foundation (ID 100010434). The fellowship code is "LCF/BQ/PR22/11920010". This study has been funded by Instituto de Salud Carlos III (ISCIII) through the project "PI21/01019" and co-funded by the European Union and by a AEI Severo Ochoa PhD fellowship (PRE2022-102586).**

## General description

This repository accompanies our [paper](https://doi.org/10.1038/s42003-025-09096-3), and provides all the tools needed to implement our proposed _Histo-μSim_ diffusion Magnetic Resonance Imaging (dMRI) technique for cancer imaging.

- [Histology to signals manual](https://github.com/radiomicsgroup/dMRIMC/blob/main/manuals/histology_to_signals.md): Guide on preparing new substrates for Monte Carlo simulations from 2D histology and generating signals
- [In silico experiment replication](https://github.com/radiomicsgroup/dMRIMC/blob/main/manuals/parameter_estimation.md): Replication of the _in silico_ experiments from our paper.

# Using Histo-μSim
Our hope is that Histo-μSim will be useful for more people in the future. However, tracing cells, processing substrates and running simulations is a both time- and resource-consuming process. Furthermore, different people will have different protocols with various bvalues and different diffusion times. To assist in the use of our tool we are also releasing some scripts along with a large signal/parameter dictionary pair, generated using a protocol with many bvalues and even more diffusion times (δ, Δ). The scripts allow you to extract a subset of the signal and parameter dictionaries, one that most closely matches the protocol your samples were scanned with, giving access to signals generated from the substrates we showcase without the need for any new simulations or signal synthesis. As this large protocol is a `PGSE`, only such protocols are supported.

## Reference protocol information
The protocol that was used to generate the signals had the following characteristics:

* `bvalues`: [0, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000]
* Gradient duration, `δ`: [5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25]
* Gradient separation, `Δ`: [5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5, 50, 52.5, 55, 57.5, 60, 62.5, 65, 70, 72.5, 80]

Note that these are the _unique_ values for each parameter. The protocol was created by making all possible combinations of the above (with `Δ` >= `δ`) leading to a total of `4761` combinations. All files describing the protocol can be found at [using_Histo_uSim/protocols/reference_protocol](https://github.com/radiomicsgroup/dMRIMC/tree/main/using_Histo_uSim/protocols/reference_protocol). Also, in an effort to escape the effects of vascular flow in the signals, which are not accounted for in our simulations, the minimum bvalue is `300` s/mm<sup>2</sup>. Both volumes that are too high and volumes that are too low are removed from the scheme and `.nii` file during processing with a tolerance of 50 s/mm2.

## Signal and parameter arrays
The signal arrays generated using this protocol are inside the [using_Histo_uSim/reference_signal_arrays](https://github.com/radiomicsgroup/dMRIMC/tree/main/using_Histo_uSim/reference_signal_arrays) folder. Similarly to what we have done for our paper, we generated the signals for 225 different cases for each substrate, 5 unique values for both intra-cellular and extra-cellular intrinsic diffusivities `D0in` and `D0ex`, and 9 for `kappa` for a total of 4050 signals. For ease of selection of cell membrane permeability, we release separate signal arrays for each permeability value, with each array having 450 different signals (5 values of `D0in`, 5 for `D0ex` times 18 substrates) as well as a numpy array with all of them combined.

Since the only variation between the signal arrays is the permeability value, all parameter arrays are the same with the exception of the value of `kappa` which is the same for all signals per array.

For this tutorial and by default, the combined signal array (all permeability values) is used, this can be changed inside the script file.

Due to file size constraints, the combined signal array could not be uploaded, to create it run `combine_arrays.py` (also creates the combined parameter array)

## Selecting array subsets
Script [get_closest_scheme.py](https://github.com/radiomicsgroup/dMRIMC/tree/main/using_Histo_uSim/get_closest_scheme.py) will find the scheme (bvalue, δ, Δ) that most closely matches the one you provide it with. In case of multiple, _consecutive_ b = 0 values created, only the first is kept and the rest are discarded with the corresponding timings also updated. The script will select the appropriate subset of the signal array you point it to wrt the protocol-describing files and save the new signal array plus a scheme file describing the protocol subset. It will also normalize the DWI niftii of your scan either via the b = 0 volume or the mean of the b = 0 volumes in case of many, and will remove the redundant b = 0 as described above.

## Selecting parameter configurations
The data we have collected span a variety of parameters describing the substrates. Running the simulations also introduces parameters such as the intrinsic diffusivities and the cell permeability. In our paper we focus on specific parameter configurations both for comparison purposes as well as the inability to just estimate all the parameters at once. When using Histo-μSim you can select which _cell size_ parameters you want to estimate via [select_parameter_configuration.py](https://github.com/radiomicsgroup/dMRIMC/tree/main/using_Histo_uSim/select_parameter_configuration.py). The available parameters are the following, for details regarding their calculation (when applicable) check the paper p.20:

* `fin`: intracellular fraction
* `mCS`: mean cell size
* `varCS`: variance of cell size
* `skewCS`: skewness of cell size
* `vCS_sph`: volume-weighted CS (vCS) index for a system with spherical geometry
* `vCS_cyl`: volume-weighted CS (vCS) index for a system with cylindrical geometry

Again, for ease of use regarding permeability, the parameter arrays are included combined and separated by the permeability value at [using_Histo_uSim/reference_param_arrays](https://github.com/radiomicsgroup/dMRIMC/tree/main/using_Histo_uSim/reference_param_arrays). By default the `select_parameter_configuration.py` script transforms the combined array.

## Example: Mice data 
As a concrete example, we show how to use Histo-μSim on the mouse data we have released on Zenodo at [https://doi.org/10.5281/zenodo.14559355](https://doi.org/10.5281/zenodo.14559355), specifically the breast cancer samples in `scans/breast`. For more information on parameter estimation check [parameter_estimation.md](https://github.com/radiomicsgroup/dMRIMC/tree/main/manuals/parameter_estimation.md) where we replicate some of the figures of our paper.

### 1. Data
To begin clone the repo and navigate to the [using_Histo_uSim](https://github.com/radiomicsgroup/dMRIMC/tree/main/using_Histo_uSim) folder:

```
git clone https://github.com/radiomicsgroup/dMRIMC.git
cd using_Histo_uSim
```

The DWI scan is already preprocessed (denoising and Gibbs unringing) so we can use it right away. We will use the following files:

* `dwi_denoise_unring_sphmean.nii`: The file containing the preprocessed and direction averaged scan
* `dwi_noise.nii`: The file containing the standard deviation of noise from MPPCA denoising
* `dwi_mask_one_sample.nii`: A mask file we have prepared covering one of specimens, a file with all of them is also included

In this case the samples were scanned in three orthogonal directions but since we are going to use the direction-averaged `.nii` file for fitting, we also have a corresponding scheme:

* `dwi_denoise_unring_sphmean.bval`
* `dwi_denoise_unring_sphmean.gdur`
* `dwi_denoise_unring_sphmean.gsep`

For the purpose of keeping this tutorial self-contained, we are including all the necessary data in the [zenodo_mouse_data](https://github.com/radiomicsgroup/dMRIMC/tree/main/using_Histo_uSim/zenodo_mouse_data) folder within the [using_Histo_uSim](https://github.com/radiomicsgroup/dMRIMC/tree/main/using_Histo_uSim) directory, where the following are already contained:

* `get_closest_scheme.py`
* `select_parameter_configuration.py`
* `mri2micro_dictml.py`
* `run_full_fitting_pipeline.py`

The `zenodo_mouse_data` should look like this:

```
zenodo_mouse_data/
├── dwi_denoise_unring_sphmean.bval
├── dwi_denoise_unring_sphmean.gdur
├── dwi_denoise_unring_sphmean.gsep
├── dwi_denoise_unring_sphmean.nii
├── dwi_mask.nii
├── dwi_mask_one_sample.nii
└── dwi_noise.nii

```


### 2. Get closest scheme, signal subset and process DWI scan file

First we need to combine the signal and parameter arrays

```
python combine_arrays.py
```

then we get the closest scheme from the available bvalues and diffusion times. We will call our target protocol `MOUSE_BREAST_EXVIVO`.


```
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
```
If all went well you should see:
```
Warning!
Max bvalue in the target scheme is 4628.46 s/mm2
Max bvalue synthesized signals are at 3000.0 s/mm2 
Volumes with bvalues higher than that are removed

-------------------------
Resulting scheme
-------------------------
[   0.  500. 2100.    0.  500. 2100.]
[11. 11. 11. 11. 11. 11.]
[16. 16. 16. 36. 36. 36.]

-------------------------
Target scheme
-------------------------
[   7.94  520.1  2063.      7.58  516.84 2056.56]
[12. 12. 12. 12. 12. 12.]
[16.5 16.5 16.5 37.  37.  37. ]

Saved closest protocol at: protocols/MOUSE_BREAST_EXVIVO

Normalizing the scan and noise, assuming that any volumes with b < 20 s/mm2 are b = 0

```
a folder `protocols` with the subfolder `MOUSE_BREAST_EXVIVO` should have been created with the files:

* `closest.bval`
* `closest.gdur`
* `closest.gsep`
* `closest.scheme`

and a file called `signal_arr_subset.npy` should have been created, containing the columns corresponding to the closest protocol.

Also, the script normalized the DWI scan and the noise by the b = 0 volume or the b = 0 mean if multiple b = 0 volumes exist. For cases when there is a small bvalue instead of 0 they are treated as b = 0 with the threshold for this controlled with the `--bval-threshold` flag. The processed DWI file should appear with the name `dwi_normalized.nii` and the noise as `dwi_noise_normalized.nii`. Finally with `--vasc-threshold 250` we set the limit for volumes with vascular signal. Here, all volumes with bvalue < 250 s/mm2 would be removed. This scan did not contain any, but in the case that it does a warning similar to the one about high bvalues will be printed and the relevant volumes will be removed.

As mentioned above you can select a subset of the available parameters, with `D0in`, `D0ex` and `kappa` always included unless a parameter array with fixed permeability is selected in which case only `D0in` and `D0ex` are included with a fixed value of `kappa`. This time we will select `fin`, `vCS_cyl`:

```
python select_parameter_configuration.py --params fin vCS_cyl
```

We now should have `param_arr_subset.npy` in the folder, with **5** columns, `fin`, `vCS_cyl`, `D0in`, `D0ex` and `kappa`

### 3. Running the fitting
Create a folder called `fitting` with

```
mkdir -v fitting
```

With all the aforementioned ingredients at hand we can now run the fitting script:

```
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
```
* `dwi_normalized.nii`: the `.nii` file we created above, normalized and with invalid volumes removed
* `signal_arr_subset.npy`: the signal array subset corresponding to the closest matching protocol
* `param_arr_subset.npy`: the parameter array subset corresponding to the closest matching protocol
* `--sldim 0`: parallelizing the fitting on the final dimension
* `--savg 3`: as mentioned above the `.nii` we are using is the average of three directions (x, y, z)
* `--ncpu 10`: using 10 threads
* `--reg "2,0.0025"`: L-norm type for the regularization, weight of the regularizer ([0.0 - 1.0])
* `--noise dwi_noise_normalized.nii`: the normalized noise file we created above
* `--mask zenodo_mouse_data/dwi_mask_one_sample.nii`: a mask file covering one of the samples
* `fitting/Histo_uSim`: the location of the output with an output string for the result files `Histo_uSim` + `_par{N}`

For more information on the script and each option check the docstring

The script will output:

```
***********************************************************************
                           mri2micro_dictml.py                         
***********************************************************************

** 4D NIFTI file with MRI measurements: dwi_normalized.nii
** MRI signal dictionary: signal_arr_subset.npy
** Tissue parameter dictionary: param_arr_subset.npy
** Algorithm for signal = f(tissue parameters) regression (forward model): rbf
** Forward regressor rbf hyperparameters: [1. 1.]
** Non-linear fitting constrained minimisation algorithm (inverse model): trust-constr
** Number of threads for parallel slice processing: 10
** Slice dimension for parallel processing: 0
** Number of words for each tissue parameter grid search: 0
** Lower bound for tissue parameters: None
** Upper bound for tissue parameters: None
** Fitting regularisation options: Lnorm = 2, weight = 0.0025
** Optional binary mask file: zenodo_mouse_data/dwi_mask_one_sample.nii
** Optional 3D NIFTI file with noise standard deviation map: dwi_noise_normalized.nii
** Number of signal averages: 3
** Output root name: fitting/Histo_uSim

    ... loading MRI measurements

    ... loading mask

    ... loading noise map

    ... training a rbf regressor to approximate the forward model signal = f(tissue parameters)
        (using 100.0 % of the dictionary)

    ... processing -- please wait

    ... saving output files
    ... done - it took 1315.9249622821808 sec

```

The `fitting` folder should contain `Histo_uSim_par{N}` with N the number of parameters selected, in the order that was specified.

The resulting map for `fin` should look like this:

<div align="center">
  <img src="https://github.com/radiomicsgroup/dMRIMC/blob/main/imgs/ex_fitting_result.png" alt="commbio" width="auto" height="auto">
</div>

### 4. Running the whole thing
We package all the above in `run_full_fitting_pipeline.sh`

```
bash run_full_fitting_pipeline.sh
```

### 5. Tips
* If you need to fix `kappa` to a specific value you need to use the appropriate pair of signal/parameter arrays from the `reference_signal_arrays`/ `reference_param_arrays` folders. The parameter arrays have `kappa` as the last column, in the case of a fixed value that column can be removed before use. You will also need to modify the `select_parameter_configuration.py` file to not automatically include the `kappa` column
* If you need to fix `D0in` or `D0ex` you can do so from the `mri2micro_dictml.py` script


## Dependencies
The above tutorial has been ran and tested with:
- **python** (3.9.23)
- **numpy** (1.26.3)
- **scipy** (1.13.1)
- **nibabel** (5.3.2)
- **matplotlib** (3.9.2)

## License
This repository is distributed under the **Attribution-NonCommercial-ShareAlike 4.0 International license** ([CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0)). Copyright (c) 2024, 2025, Fundació Privada Institut d’Investigació Oncològica de Vall d’Hebron (Vall d'Hebron Institute of Oncology (VHIO), Barcelona, Spain). All rights reserved. Link to license [here](https://github.com/radiomicsgroup/dMRIMC/blob/main/license.txt). 
