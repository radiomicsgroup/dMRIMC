## Simulation-informed microstructure parameter estimation
Signals generated from histology-derived substrates can be used to inform microstructural parameter estimation, by developing numerical signal model that output diffusion MRI (dMRI) signals given tissue parameters, for any dMRI protocol of interest. These numerical signal models can be embedded into standard maximum likelihood fitting, improving the performances compared to fitting standard analytical signal expressions. This is due to the fact that the latter often reliy on naive modelling idealisations of the biology of the tissues being imaged (e.g., regular cylinders/spheres for cells). 

This tutorial provides an example of simulation-informed parameter estimation on signals generated from the histologically realistic cancer substrates included in our [preprint](https://www.medrxiv.org/content/10.1101/2024.07.15.24310280v1). In doing so, we reproduce some of the figures reporing on the _in silico_ experiments of our paper.

**Under development! Everything will be here in the coming days -- apologies for the inconvenience!**

## Data
We provide you already with the synthetic signals and corresponding tissue parameters generated from our 18 histology-derived substrates. These have been stored in NIFTI format, following the same steps described in the previous tutorial [here](https://github.com/radiomicsgroup/dMRIMC/blob/main/manuals/histology_to_signals.md). For each substrate, we use 10 unique values for both intra-cellular and extra-cellular intrisnic diffusivities $`D_{0|in}$ and $D_{0|ex}$, for a total of 1800 signals. The substrates are shown here:

FIG

We have partitioned them according to a leave-one-out procedure: we use 17 out of 18 substrates to learn a numerical forward model, which then we fit on the 18th substrate. The forward model is learnt via radial basis function regression of noise-free signals, while fitting (i.e., model inversion) is performed on signals corrupted with Rician noise (SNR = 50 and 20). Leave-one-out partitioned signals and tissue parameters are stored respectively in the `leave_one_out/signal_np_files` and `leave_one_out/param_np_files` folders.

Fitting accounts for Rician bias, and a noise estimate is obtained with `dwidenoise` from [MRtrix3](https://mrtrix.readthedocs.io/en/latest/reference/commands/dwidenoise.html). Fitted parameters are also stored in NIFTI format in `leave_one_out/niftiis` (MC-informed fitting signals) and `leave_one_out/analytical` (fitting of a standard analytical signal model analytical expression). 

Here we consider three different acquisition protocols, which are the same we used in our preprint. These are: 
* `PGSEin`: a pulsed-gradient spin echo (PGSE) protocol, which in the paper we used for _in vivo_ imaging (salient characteristics: XXX)
* `PGSEex`: a second PGSE protocol, which in the paper we used for _ex vivo_ imaging;
* `TRSE`: a diffusion-weighted (DW) twice-refocussed spin echo (TRSe) protocol, which  in the paper we also used for _in vivo_ imaging.


## Fitting configurations
We fit two different simulation-informed signal models. Due to practical code implementation, we refer to these as two different _fitting configurations_. These correpsond to forward models 1 and forward 2 in our paper, namely:

- in forward model 1 (fitting configuration number 9 in our code) we estimate $`f_{in},`$ vCS, $`D_{0|in}, D_{0|ex}`$;

- in forward model 2 (fitting configuration number 13 in our code) we estimate $`f_{in},`$ mCS, varCS, skewCS, $`D_{0|in}, D_{0|ex}`$.

The parameters in forward model 1 were chosen so that we could compare the performance of simulation-informed fitting vs a well-established two-compartment model, accounting for restricted diffusion within cylinders and hindered extra-cellular Gaussian diffusion. Remember that our substrates where obtained on 2D histology, and thus feature cylndrical symmetry (that is why we used cylinders, instead of spheres).

Forward model 2 instead investigates the feasibility of data-driven cell size distribution characterisation, where its first three moments are estimated without imposing any parametric form (e.g., the Gamma distribution).

## Fitting procedures
1. Fitting the analytical model (check script help for more information). This script performs the fitting of the analytical expression using all the signals at the same time.
    - `python dri2mc_maxlikcyl.py --noise analytical/noise_maps/PGSEin_noise_analytical_SNR50_all_signals.nii --modstr DinDex --pmin 8.0,0.8,0.0,0.5 --pmax 20.0,3.0,0.9,3.0 --sldim 0 --nw 12 --ncpu 10 analytical/niftiis/PGSEin_all_signals_config_9_SNR_50.nii PGSEin_for_maxlik.bval res_maxlik/old_way_cylinders/config_9/SNR_50`

2. Fitting using the MC-generated signals with the PGSEin protocol. The `run_all.py` script will run the fitting routine for all protocols and configurations.
    - `python run_all.py`

To reproduce Fig 5 from our paper i.e contour plots comparing the MC-informed parameter estimation with the analyical expression for a hindered-restricted diffusion model use `explore_niftii_MC_informed.py` and `explore_niftii_ANALYTICAL.py` after the fittings.

1. `python explore_niftii_MC_informed.py PGSEin 50 9`: This will create Fig 5a (MC-informed fitting)
2. `python explore_niftii_ANALYTICAL.py 50`: This will create Fig 5b (analytical fitting)


## Code
### Fitting
 - `dri2mc_maxlikcyl.py`: Script for fitting the analytical model. Requires a noise map of the signal array, an array containing all signals corrupted with Rician noise and a text file with the bvalues of the protocol. All files but the bvalue specifications must be in niftii format. 
   
 - `mri2micro_dictml.py`: Script for the fitting using the MC-generated signals (an updated version can be found at [bodymritools](https://github.com/fragrussu/bodymritools)). Requires a noise map of the signal array, a signal array corrupted with Rician noise generated from a substrate and two numpy arrays, one with the signals that will be used for training and one containing the correspoding parameters of the substrate.
   
 - `run_all.py`: Convenience script that runs the MC-informed fitting for all cases (all protocols and all configurations/forward models)

### Plotting
- `explore_niftii_MC_informed.py`: This script will generate figures based on the selected protocol, SNR and configuration. The figures are contour density plots that show the quality of the fitting compared to the ground truth.

<div align="center">
    <img src="https://github.com/radiomicsgroup/dMRIMC/blob/main/parameter_estimation/figs/PGSEin_SNR_50_config_9_all_params.jpg" alt="mcinformed" width="auto" height="auto">
</div>

- `explore_niftii_ANALYTICAL.py`: Similar to the one above but only takes SNR as an argument. Generates the plot of the analytical fitting compared to the ground truth

<div align="center">
    <img src="https://github.com/radiomicsgroup/dMRIMC/blob/main/parameter_estimation/figs/PGSEin_old_way_cylinders_SNR_50_config_9_all_params_ANALYTICAL.jpg" alt="mcinformed" width="auto" height="auto">
</div>
