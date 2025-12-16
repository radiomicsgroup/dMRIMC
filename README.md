## Histology-informed microstructural diffusion simulations for MRI cancer characterisation — the Histo-μSim framework

<div align="center">
  <img src="https://github.com/radiomicsgroup/dMRIMC/blob/main/imgs/commbio25.png" alt="commbio" width="auto" height="auto">
</div>

If you use the code released in this repository in your research, please cite our article:

Athanasios Grigoriou, et al. **"Histology-informed microstructural diffusion simulations for MRI cancer characterisation — the Histo-μSim framework"**. Communications Biology 2025, 8: 1695, doi: [10.1101/2024.07.15.24310280](https://doi.org/10.1038/s42003-025-09096-3).

This repository was developed by Athanasios Grigoriou (<agrigoriou@vhio.net>) and Francesco Grussu (<fgrussu@vhio.net>). **The project that gave rise to these results received the support of a fellowship from ”la Caixa” Foundation (ID 100010434). The fellowship code is "LCF/BQ/PR22/11920010". This study has been funded by Instituto de Salud Carlos III (ISCIII) through the project "PI21/01019" and co-funded by the European Union and by a AEI Severo Ochoa PhD fellowship (PRE2022-102586).**

## General description

This repository accompanies our [paper](https://doi.org/10.1038/s42003-025-09096-3), and provides all the tools needed to implement our proposed _Histo-μSim_ diffusion Magnetic Resonance Imaging (dMRI) technique for cancer imaging. 

Below we show you **how to start using Histo-μSim immediately on your own data**. Additionally, we have also written two more tutorials:

- [Histology-to-signals manual](https://github.com/radiomicsgroup/dMRIMC/blob/main/manuals/histology_to_signals.md): to guide you on preparing new substrates for Monte Carlo simulations from 2D histology, in case you want to run more simulations;
- [In silico experiment replication](https://github.com/radiomicsgroup/dMRIMC/blob/main/manuals/parameter_estimation.md): to show you how to replicate some of the _in silico_ experiments performed in our paper.

# Using Histo-μSim
Our hope is that Histo-μSim will be useful for more people in the future. However, we are aware that tracing cells, processing substrates and running simulations is a time- and resource-consuming process, and that it may result complicated for clinical labs focussing on applied imaging. 

To assist in the use of our tool, we have prepared a **rich dictionary of synthetic signals that you can download and deploy immediately to fit Histo-μSim on your own diffusion MRI scans**. 

This signal dictionary corresponds to a very rich protocol with multiple b-values and even more diffusion times (δ, Δ), within which you will certainly find the protocol that you used to acquire your own data. The dictionary comes with a set of scripts that allow you to extract the subset of the synthetic signals that most closely matches the protocol that you have acquired. This will give you access to the full potential of Histo-μSim, without the need for any new simulations or signal synthesis. 


## Rich protocol information
We have generated synthetic signals for a very rich protocol where you will be able to find your own diffusion measurements with a precision of as few as 50 s/mm<sup>2</sup> for $b$, and just 2.5 ms for δ and Δ. The synthetic signals were obtained for all possible combinations of

* $b$ = [0, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000] s/mm<sup>2</sup>
* Gradient duration δ = [5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25] ms
* Gradient separation Δ = [5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5, 50, 52.5, 55, 57.5, 60, 62.5, 65, 70, 72.5, 80] ms

making sure of course that Δ $\geq$ δ. This leads to a total of `4761` combinations with unique ($b$,δ,Δ). All files describing the protocol can be found at [using_Histo_uSim/protocols/reference_protocol](https://github.com/radiomicsgroup/dMRIMC/tree/main/using_Histo_uSim/protocols/reference_protocol). Importantly, note that the minimum bvalue is `300` s/mm<sup>2</sup>: this choice allows us to minimise the contribution of the vascular signal _in vivo_, since our simulations do not account for capillary perfusion. Also, our code will exclude from the fitting diffusion measurements that you might have acquired for b-values that are higher/lower than the maximum/minimum b-value that we have simulated. 


Note that **only pulsed-gradient spin echo (PGSE) protocols are supported**. The synthetic signals stored in the [reference_protocol](https://github.com/radiomicsgroup/dMRIMC/tree/main/using_Histo_uSim/protocols/reference_protocol) folder are in the format of NumPy matrices where rows represents different microstructures, while columns different measurements in the dMRI protocol.


## Signal and parameter arrays
The signal arrays generated using this protocol are inside the [using_Histo_uSim/reference_signal_arrays](https://github.com/radiomicsgroup/dMRIMC/tree/main/using_Histo_uSim/reference_signal_arrays) folder. **We have generated signals for the 18 cancer substrates that we produced for our paper**. The substrates can be downloaded from Grigoriou et al, Zenodo 2024, [doi: 10.5281/zenodo.14559103](https://doi.org/10.5281/zenodo.14559103).

We have generated signals for **225 unique realisations of each substrate**, obtained by varying the intrinsic intra-cellular diffusivity `D0in` (5 values), the intrinsic extra-cellular diffusivity `D0ex` (5 values), and the cell membrane permeability `kappa` (9 values). This leads to a **total of 4050 signals for each ($b$,δ,Δ) measurement**. Note that we release separate signal/parameter arrays for each `kappa` value, with each array having 450 different signals (5 values of `D0in`, 5 for `D0ex` times 18 substrates) due to file size constraints (see the [reference_signal_arrays](https://github.com/radiomicsgroup/dMRIMC/tree/main/using_Histo_uSim/reference_signal_arrays) and the [reference_param_arrays](https://github.com/radiomicsgroup/dMRIMC/tree/main/using_Histo_uSim/reference_param_arrays) folders).

The range of variation for `D0in`, `D0ex` and `kappa` are:
* `D0in`: from 0.8 μm<sup>2</sup>/ms to 3.0 μm<sup>2</sup>/ms
* `D0ex`: from 0.8 μm<sup>2</sup>/ms to 3.0 μm<sup>2</sup>/ms
* `kappa`: from 0 μm/s to 40 μm/s

## Selecting signal subsets and tissue parameter configurations
Below you will find a practical example that will illustrate how to use these synthetic signals to fit _Histo-μSim_ on your data. 

Briefly, script [get_closest_scheme.py](https://github.com/radiomicsgroup/dMRIMC/tree/main/using_Histo_uSim/get_closest_scheme.py) will find subset of synthetic measurements that most closely match the acquisition scheme ($b$, δ, Δ) that you used to acquire your data. Additionally, it will also normalize the dMRI measurements you provided in the input NIFTI file so that the signal at $b$ = 0 is 1 (our synthetic signals are bound within 0 and 1).

Once the subset of synthetic measurements is found and your dMRI scan has been normalised, you will have to choose which tissue parameters to fit using the [select_parameter_configuration.py](https://github.com/radiomicsgroup/dMRIMC/tree/main/using_Histo_uSim/select_parameter_configuration.py). The available parameters are the following (for details regarding their calculation page 20 of our paper):

* `fin`: intracellular fraction (normalised), ranging from 0.023 to 0.867
* `mCS`: mean cell size in μm, ranging from 6.1 μm to 15.9 μm
* `varCS`: variance of cell size in μm<sup>2</sup>, ranging from 2.3 μm<sup>2</sup> to 19.7 μm<sup>2</sup>
* `skewCS`: skewness of cell size (dimensionless), ranging from -0.52 to 0.86
* `vCS_sph`: volume-weighted CS (vCS) index for a system with spherical geometry in μm, ranging from 8.2 μm to 19.9 μm 
* `vCS_cyl`: volume-weighted CS (vCS) index for a system with cylindrical geometry in μm, ranging from 7.9 μm to 19.1 μm

Tissue parameters corresponding to synthetic signals are stored in the folder  [using_Histo_uSim/reference_param_arrays](https://github.com/radiomicsgroup/dMRIMC/tree/main/using_Histo_uSim/reference_param_arrays) as NumPy matrices where rows represent different microstructure realisations, while columns are the tissue parameters corresponding to each microstructure realisation.  

## Example: mouse data 
As a concrete example, we show how to use Histo-μSim on the mouse data we have released on Zenodo at [https://doi.org/10.5281/zenodo.14559355](https://doi.org/10.5281/zenodo.14559355), specifically the breast cancer samples in `scans/breast`. For more information on parameter estimation check [parameter_estimation.md](https://github.com/radiomicsgroup/dMRIMC/tree/main/manuals/parameter_estimation.md) where we replicate some of the figures of our paper.

### 1. Data
To begin, clone the Histo-μSim repo and navigate to the [using_Histo_uSim](https://github.com/radiomicsgroup/dMRIMC/tree/main/using_Histo_uSim) folder:

```
git clone https://github.com/radiomicsgroup/dMRIMC.git
cd dMRIMC
cd using_Histo_uSim
```

The folder [zenodo_mouse_data](https://github.com/radiomicsgroup/dMRIMC/tree/main/using_Histo_uSim/zenodo_mouse_data) within [using_Histo_uSim](https://github.com/radiomicsgroup/dMRIMC/tree/main/using_Histo_uSim) already contains a pre-processed dMRI scan, ready to use for this tutorial (so there is no need for you to download anything). Pre-processing included MP-PCA denoising, Gibbs unringing and directional averaging. We will use the following files:

* `dwi_denoise_unring_sphmean.nii`: The file containing a preprocessed and directionally averaged scan
* `dwi_noise.nii`: The file containing the standard deviation of noise from MPPCA denoising
* `dwi_mask_one_sample.nii`: A mask file we have prepared covering one of specimens (a file with all of them is also included)

The scan comes from our [Zenodo](https://doi.org/10.5281/zenodo.14559355) data set, which contains the mouse data used in our Histo-μSim paper. The data was acquired in tissue samples were scanned acquiring diffusion images along three orthogonal directions. However, note that **Histo-μSim is meant to be used on directionally-averaged signals** - that is the reason why we will be studying the 4D `dwi_denoise_unring_sphmean.nii` file. The b-values and gradient timing corresponding to each volume of `dwi_denoise_unring_sphmean.nii` are indicated in the following additional files: 

* `dwi_denoise_unring_sphmean.bval`(b-values in s/mm<sup>2</sup>)
* `dwi_denoise_unring_sphmean.gdur` (gradient duration δ in ms)
* `dwi_denoise_unring_sphmean.gsep` (gradient separation Δ in ms)

In practice, the `zenodo_mouse_data` should look like this:

```
zenodo_mouse_data/
├── dwi_denoise_unring_sphmean.bval
├── dwi_denoise_unring_sphmean.gdur
├── dwi_denoise_unring_sphmean.gsep
├── dwi_denoise_unring_sphmean.nii
├── dwi_denoise_unring_sphmean_difftimefixed.bval
├── dwi_denoise_unring_sphmean_difftimefixed.gdur
├── dwi_denoise_unring_sphmean_difftimefixed.gsep
├── dwi_denoise_unring_sphmean_difftimefixed.nii
├── dwi_mask.nii
├── dwi_mask_one_sample.nii
└── dwi_noise.nii

```

The folder [using_Histo_uSim](https://github.com/radiomicsgroup/dMRIMC/tree/main/using_Histo_uSim) also contains all the necessary code: 

* `get_closest_scheme.py`
* `select_parameter_configuration.py`
* `mri2micro_dictml.py`
* `run_full_fitting_pipeline.py`



### 2. Get synthetic signal subsets and pre-process the input dMRI scan

First we need to combine the signal and parameter arrays. This step is necessary because we uploaded our signal/parameter dictionaries through multiple files to comply with GitHub restrictions of file size, and is done simply by running 

```
python combine_arrays.py
```


Afterwards, we need to extract the closest scheme from the available b-values and diffusion times using [`get_closest_scheme.py`](https://github.com/radiomicsgroup/dMRIMC/blob/main/using_Histo_uSim/get_closest_scheme.py). We will call our target protocol `MOUSE_BREAST_EXVIVO` (but you can use any name you like):


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
* `dwi_noise_normalized.nii`
* `dwi_normalized.nii`
* `signal_arr_subset.npy`

The `signal_arr_subset.npy` and  file should have been created, containing the columns corresponding to the closest protocol. Also, the script will have normalized the input dMRI scan `zenodo_mouse_data/dwi_denoise_unring_sphmean.nii` scan and the noise level `zenodo_mouse_data/dwi_noise.nii` by the mean b = 0 volume. Do not worry if you do not have a noise map, as it is an optional input parameter. For cases when there is a small bvalue instead of 0 they are treated as b = 0 with the threshold for this controlled with the `--bval-threshold` flag. The processed dMRI file should appear with the name `dwi_normalized.nii` and the noise as `dwi_noise_normalized.nii`. 

Note that in the example above we used the option `--vasc-threshold 250`. This sets a threshold to discard b-values lower than the threshold, and is meant to remove dMRI measurements with non-negligible vascular signal contributions. In this case, all volumes with bvalue < 250 s/mm<sup>2</sup> (except for the b = 0 volume) would be removed. This scan did not contain any, but in the case that it does a warning similar to the one about high bvalues will be printed and the relevant volumes will be removed.

### 3. Select which tissue parameters to estimate

As mentioned above, you can now select a subset of the available parameters to fit using [`select_parameter_configuration.py`](https://github.com/radiomicsgroup/dMRIMC/blob/main/using_Histo_uSim/select_parameter_configuration.py). Note that **`D0in`, `D0ex` and `kappa` are always included**. We will focus on the estimation of `fin` and `vCS_cyl`, exactly as we did in our paper:

```
python select_parameter_configuration.py --params fin vCS_cyl --output-folder MOUSE_BREAST_EXVIVO
```

You will see that this will have created the tissue parameter file `param_arr_subset.npy` (**5** columns: `fin`, `vCS_cyl`, `D0in`, `D0ex` and `kappa`), corresponding with the synthetic signals array `signal_arr_subset.npy`. 

### 4. Perform Histo-μSim model fitting
We are now ready to use the synthetic signals and the corresponding tissue parameter file for Histo-μSim fitting. Let's create a folder called `fitting` to store the fitting results:

```
mkdir -v protocols/MOUSE_BREAST_EXVIVO/fitting
```

With all the aforementioned ingredients at hand we can now run the fitting script using the [`mri2micro_dictml.py`](https://github.com/radiomicsgroup/dMRIMC/blob/main/using_Histo_uSim/mri2micro_dictml.py) script:

```
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
```

We are distributing a copy of `mri2micro_dictml.py` here, but note that it was originally released as part of the [BodyMRItools](https://github.com/fragrussu/bodymritools/) repository. We invite you to check it out, as it include many more scripts that can be of help to work with body diffusion imaging!

A quick comment on the inputs taken by `mri2micro_dictml.py`:
* `dwi_normalized.nii`: the `.nii` file we created above, normalized and with measurements that fall outside the range of the simulated protocol removed
* `signal_arr_subset.npy`: the signal dictionary corresponding to the closest matching protocol
* `param_arr_subset.npy`: the tissue parameter dictionary corresponding to the closest matching protocol
* `--sldim 0`: spatial dimension along which paralleling the fitting (0: image dimension i; 1: image dimension j; 2: image dimension k)
* `--savg 3`: number of signal averages used to acquire the input `.nii` above. We are using 3 as the scan is the average of three directions (x, y, z)
* `--ncpu 10`: number of threads to be used for parallel computing (10 threads in this example)
* `--reg "2,0.0025"`: use regularisation for model witting. We use an L2 regularisation with regularisation weight of 0.0025 (choose a number betwee 0, for no regularisation, and 1) 
* `--noise dwi_noise_normalized.nii`: the normalized noise file we created above, to model noise floor bias
* `--mask zenodo_mouse_data/dwi_mask_one_sample.nii`: a mask file covering one of the samples
* `fitting/Histo_uSim`: the location of the output with an output string for the result files `Histo_uSim` + `_par{N}.nii`

For more information on all the options of `mri2micro_dictml.py`, simply type `python python mri2micro_dictml.py -h`. Running the script will output:

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

The `fitting` folder should contain `Histo_uSim_par{N}` with N the number of parameters selected, in the order that was specified: parameter 1 is `fin`, parameter 2 is `vCS_cyl` (in μm), parameter 3 is `D0in` (in μm<sup>2</sup>/ms), parameter 4 is `D0ex` (in μm<sup>2</sup>/ms), parameter 5 is `kappa` (in μm/s). 

The resulting map for `fin` should look like this:

<div align="center">
  <img src="https://github.com/radiomicsgroup/dMRIMC/blob/main/imgs/ex_fitting_result.png" alt="commbio" width="auto" height="auto">
</div>

### 5. A script to run the whole thing at once
We have packaged all the above in one command line script, called [`run_full_fitting_pipeline.sh`](https://github.com/radiomicsgroup/dMRIMC/blob/main/using_Histo_uSim/run_full_fitting_pipeline.sh). You can run it by typing:

```
bash run_full_fitting_pipeline.sh
```

### 6. Additional tips
All the .py python scripts described in this tutorial have their own help manual. To check it, simply type
```
python <SCRIPT>.py -h
```
(for example, `python mri2micro_dictml.py -h`).


The default implementation of Histo-μSim attempts to resolve diffusion restriction lengths and cell permeability at once, and hence requires dMRI protocols featuring multiple diffusion times. **If your protocol contains b-values acquired using only a single diffusion time, mathematically you cannot resolve all these properties together, and you might need to fix some of the parameters to specific _and hoc_ values**. 


However, be aware: **fixing tissue parameters will bias the estimation of the other free parameters!** You would then need to be extremely careful when interpreting the output maps.


The first thing you can do is to fix `kappa` to a specific value. Additionally, you would need to fix  `D0in` and, potentially, `D0ex`, depending on the number of b-values that you have acquired. The script [`run_full_fitting_pipeline_fixpars.sh`](https://github.com/radiomicsgroup/dMRIMC/blob/main/using_Histo_uSim/run_full_fitting_pipeline_fixpars.sh) shows you how to do it in practice.



## Dependencies
The above tutorial has been ran and tested with:
- **python** (3.9.23)
- **numpy** (1.26.3)
- **scipy** (1.13.1)
- **nibabel** (5.3.2)
- **matplotlib** (3.9.2)

## License
This repository is distributed under the **Attribution-NonCommercial-ShareAlike 4.0 International license** ([CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0)). Copyright (c) 2024, 2025, Fundació Privada Institut d’Investigació Oncològica de Vall d’Hebron (Vall d'Hebron Institute of Oncology (VHIO), Barcelona, Spain). All rights reserved. Link to license [here](https://github.com/radiomicsgroup/dMRIMC/blob/main/license.txt). 
