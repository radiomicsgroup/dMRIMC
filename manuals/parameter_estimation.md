# Simulation-informed microstructure parameter estimation
Signals generated from histology-derived substrates can be used to inform microstructural parameter estimation, by developing numerical signal model that output diffusion MRI (dMRI) signals given tissue parameters, for any dMRI protocol of interest. These numerical signal models can be embedded into standard maximum likelihood fitting, improving the performances compared to fitting standard analytical signal expressions. This is due to the fact that the latter often reliy on naive modelling idealisations of the biology of the tissues being imaged (e.g., regular cylinders/spheres for cells). 

This tutorial provides an example of simulation-informed parameter estimation on signals generated from the histologically realistic cancer substrates included in our [preprint](https://www.medrxiv.org/content/10.1101/2024.07.15.24310280v1). In doing so, we reproduce some of the figures reporing on the _in silico_ experiments of our paper.

**Under development! Everything will be here in the coming days -- apologies for the inconvenience!**

## Data
We provide you already with the synthetic signals and corresponding tissue parameters generated from our 18 histology-derived substrates. These have been stored in NIFTI format, following the same steps described in the previous tutorial [here](https://github.com/radiomicsgroup/dMRIMC/blob/main/manuals/histology_to_signals.md). For each substrate, we use 10 unique values for both intra-cellular and extra-cellular intrisnic diffusivities $D_{0|in}$ and $D_{0|ex}$, for a total of 1800 signals. The substrates are shown here:

<p align="center">
  <img src="https://github.com/radiomicsgroup/dMRIMC/blob/main/imgs/all_subs_with_metrics.jpg" width="80%" />
</p>

We have partitioned them according to a leave-one-out procedure: we use 17 out of 18 substrates to learn a numerical forward model, which then we fit on the 18th substrate. The forward model is learnt via radial basis function regression of noise-free signals, while fitting (i.e., model inversion) is performed on signals corrupted with Rician noise (SNR = 50 and 20). Leave-one-out partitioned signals and tissue parameters are stored respectively in the `leave_one_out/signal_np_files` and `leave_one_out/param_np_files` folders.

Fitting accounts for Rician bias, and a noise estimate is obtained with `dwidenoise` from [MRtrix3](https://mrtrix.readthedocs.io/en/latest/reference/commands/dwidenoise.html). Fitted parameters are also stored in NIFTI format in `leave_one_out/niftiis` (MC-informed fitting signals) and `leave_one_out/analytical` (fitting of a standard analytical signal model analytical expression). 

Here we consider three different acquisition protocols, which are the same we used in our preprint. These are: 
* `PGSEin`: a pulsed-gradient spin echo (PGSE) protocol, which in the paper we used for _in vivo_ imaging

(salient characteristics: 3 b = 0 and 18 DW measurements: b = {50, 100, 400, 900, 1200, 1500, 50, 100, 400, 900, 1200, 1500, 50, 100, 400, 900, 1200, 1500} s/mm2 , δ = {3.9, 5.2, 9.2, 15.0, 18.2, 21.0, 3.9, 5.2, 9.2, 13.0, 15.8, 18.5, 3.9, 5.2, 9.2, 13.0, 15.8, 18.5} ms, Δ = {27.8, 29.0, 33.0, 28.7, 31.8, 34.7, 7.8, 29.0, 33.0, 37.0, 39.6, 42.3, 7.8, 29.0, 33.0, 37.0, 39.6, 42.3} ms)

* `PGSEex`: a second PGSE protocol, which in the paper we used for _ex vivo_ imaging 

(salient characteristics: 2 b = 0 and 6 DW measurements: b = {0, 500, 2000, 4500} s/mm2 acquired for each of Δ = {16.5, 37.0} ms, with δ = 12 ms)

* `TRSE`: a diffusion-weighted (DW) twice-refocussed spin echo (TRSe) protocol, which  in the paper we also used for _in vivo_ imaging 

(salient characteristics: 3 b = 0 and 18 DW measurements: b = {0, 50, 100, 400, 900, 1200,
1600} s/mm2 , repeated for 3 different diffusion times. The duration/separation of the gradient lobes for the 3 diffusion times were: δ1 = {8.9, 13.2, 18.9} ms, δ2 = {17.6, 19.3, 21.0} ms, δ3 = {20.4, 24.8, 30.5} ms, δ4 = {6.0, 7.7, 9.5} ms, Δ1,2 = {17.4, 21.7, 27.5} ms, Δ1,4 = {63.9, 74.2, 87.5} ms)


## Fitting configurations
We fit two different simulation-informed signal models. Due to practical code implementation, we refer to these as two different _fitting configurations_. These correpsond to forward models 1 and forward 2 in our paper, namely:

- in forward model 1 (fitting configuration number 9 in our code) we estimate $`f_{in},`$ vCS, $`D_{0|in}, D_{0|ex}`$;

- in forward model 2 (fitting configuration number 13 in our code) we estimate $`f_{in},`$ mCS, varCS, skewCS, $`D_{0|in}, D_{0|ex}`$.

Above, the parameters have the following meaning:
* $f_{in}$: intra-cellular signal fraction
* $vCS$: characteristic volume-weighted cell diameter
* $mCS$: mean cell diameter
* $varCS$: cell diameter variance
* $skewCS$: cell diameter skewness
* $D_{0|in}$: intrinsic intra-cellular diffusivity
* $D_{0|ex}$: intrinsic extra-cellular diffusivity.

The parameters in forward model 1 were chosen so that we could compare the performance of simulation-informed fitting vs a well-established two-compartment model, accounting for restricted diffusion within cylinders and hindered extra-cellular Gaussian diffusion. Remember that our substrates where obtained on 2D histology, and thus feature cylndrical symmetry (that is why we used cylinders, instead of spheres). Note that the analytical model can only be fitted on PGSE data, and not on TRSE.

Forward model 2 instead investigates the feasibility of data-driven cell size distribution characterisation, where its first three moments are estimated without imposing any parametric form (e.g., the Gamma distribution).


### Monte Carlo simulation-informed fitting
We provide you with the code to perform model fitting on the simulated signals. 

The script `run_all.py`: this script performs the MC-informed fitting for all cases (all protocols, i.e., `PGSEin`, `PGSEex` and `TRSE`; and both forward models 1 and 2). This script relies on the `mri2micro_dictml.py` routine (note that this is a slightly older version compared to [mri2micro_dictml.py](https://github.com/fragrussu/bodymritools/blob/main/mrifittools/mri2micro_dictml.py), released in the [bodymritools](https://github.com/fragrussu/bodymritools) python repository). To run it, simply navigate to `parameter_estimation/leave_one_out/` and run the script:

```
cd parameter_estimation/leave_one_out/
python run_all.py
```


### Fitting an analytical signal model
Conversely, you can use the `dri2mc_maxlikcyl.py` script for fitting forward model 1 on the PGSE protocols. The script requires the noisy signals to fit in NIFTI format, a text file with the diffusion protocol and an optional noise map, also in NIFTI format. For example, you can fit this analytical model on signals generated according to the `PGSEin` protocol like this:

```
cd parameter_estimation/leave_one_out/
python dri2mc_maxlikcyl.py --noise analytical/noise_maps/PGSEin_noise_analytical_SNR50_all_signals.nii --modstr DinDex --pmin 8.0,0.8,0.0,0.5 --pmax 20.0,3.0,0.9,3.0 --sldim 0 --nw 12 --ncpu 10 analytical/niftiis/PGSEin_all_signals_config_9_SNR_50.nii PGSEin_for_maxlik.bval res_maxlik/old_way_cylinders/config_9/SNR_50
```

Fitting this analytical two-pool model provides estimates of:
* $f_{in}$: intra-cellular signal fraction
* $vCS$: characteristic volume-weighted cell diameter
* $ADC_{ex}$: extra-cellular apparent diffusion coefficient.
   


### Plotting fitting results
We also provide you with scripts to generate scatter density plots that correlate estimated vs ground truth tissue parameters. We include two scripts:

- `explore_niftii_MC_informed.py`: results from MC-informed fitting for a desired protocol, SNR and forward model. You can use it like this:
```
cd parameter_estimation
python explore_niftii_MC_informed.py PGSEin 50 9
python explore_niftii_MC_informed.py TRSE 20 13
python explore_niftii_MC_informed.py PGSEex 50 13

```

- `explore_niftii_ANALYTICAL.py`: results from fitting the analytical signal model protocol `PGSEin`. It requires as input the SNR for which results should be shown, i.e., 
```
cd parameter_estimation
python explore_niftii_ANALYTICAL.py 50
python explore_niftii_ANALYTICAL.py 20
```

With these scripts, you can generate the plots related to the _in silico_ parameter estimation experiments of our paper. For example, results from fitting foward model 1 on the `PGSEin` protocol with MC-informed fitting ... 

<div align="center">
    <img src="https://github.com/radiomicsgroup/dMRIMC/blob/main/parameter_estimation/figs/PGSEin_SNR_50_config_9_all_params.jpg" alt="mcinformed" width="auto" height="auto">
</div>

... or with a two-comparment analytical signal model: 

<div align="center">
    <img src="https://github.com/radiomicsgroup/dMRIMC/blob/main/parameter_estimation/figs/PGSEin_old_way_cylinders_SNR_50_config_9_all_params_ANALYTICAL.jpg" alt="mcinformed" width="auto" height="auto">
</div>

