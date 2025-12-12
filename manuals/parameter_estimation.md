# Simulation-informed microstructural parameter estimation
Signals generated from histology-derived substrates can be used to inform microstructural parameter estimation, by developing numerical signal models that output diffusion MRI (dMRI) signals given tissue parameters, for any dMRI protocol of interest. These numerical signal models can be embedded into standard maximum likelihood fitting, improving the performances compared to fitting standard analytical signal expressions. This is due to the fact that the latter often rely on modelling idealisations of the biology of the tissues being imaged (e.g., regular cylinders/spheres for cells). 

This tutorial provides an example of simulation-informed parameter estimation on signals generated from the histologically realistic cancer substrates included in our [paper](https://doi.org/10.1038/s42003-025-09096-3) and released on [Zenodo](https://doi.org/10.5281/zenodo.14559103). In doing so, we reproduce some of the figures reporing on the _in silico_ experiments of our paper.


## Data
We provide you already with the synthetic signals and corresponding tissue parameters generated from our 18 histology-derived substrates. These have been stored in NIFTI and npy format, following the same steps described in the previous tutorial [here](https://github.com/radiomicsgroup/dMRIMC/blob/main/manuals/histology_to_signals.md). For each substrate, we use 5 unique values for both intra-cellular and extra-cellular intrinsic diffusivities D<sub>0|in</sub> and D<sub>0|ex</sub>, and 9 for $\kappa$ for a total of 4050 signals. The substrates are shown here:

<p align="center">
  <img src="https://github.com/radiomicsgroup/dMRIMC/blob/main/imgs/all_subs_with_metrics.jpg" width="80%" />
</p>

We have partitioned them according to a leave-one-out procedure: we use 17 out of 18 substrates to learn a numerical forward model, which then we fit on the 18th substrate. The forward model is learnt via radial basis function regression of noise-free signals, while fitting (i.e., model inversion) is performed on signals corrupted with Rician noise (SNR = 50 and 20). Leave-one-out partitioned signals and tissue parameters are stored respectively in the [`parameter_estimation/LeaveOneOut`](https://github.com/radiomicsgroup/dMRIMC/tree/main/parameter_estimation/LeaveOneOut) folder.

Fitting accounts for Rician bias, and a noise estimate is obtained with `dwidenoise` from [MRtrix3](https://mrtrix.readthedocs.io/en/latest/reference/commands/dwidenoise.html). Fitted parameters are also stored in NIFTI format in [`parameter_estimation/dict_fit`](https://github.com/radiomicsgroup/dMRIMC/tree/main/parameter_estimation/dict_fit) (MC-informed fitting signals).

In the paper we consider three different acquisition protocols, here we describe one of them, `PGSEin`, as an example:
* `PGSEin`: a pulsed-gradient spin echo (PGSE) protocol, which in the paper we used for _in vivo_ imaging

(salient characteristics: 3 b = 0 and 18 DW measurements: b = {50, 100, 400, 900, 1200, 1500, 50, 100, 400, 900, 1200, 1500, 50, 100, 400, 900, 1200, 1500} s/mm<sup>2</sup> , δ = {3.9, 5.2, 9.2, 15.0, 18.2, 21.0, 3.9, 5.2, 9.2, 13.0, 15.8, 18.5, 3.9, 5.2, 9.2, 13.0, 15.8, 18.5} ms, Δ = {27.8, 29.0, 33.0, 28.7, 31.8, 34.7, 7.8, 29.0, 33.0, 37.0, 39.6, 42.3, 7.8, 29.0, 33.0, 37.0, 39.6, 42.3} ms)

This figure illustrates a PGSE sequence:

<p align="center">
  <img src="https://github.com/radiomicsgroup/dMRIMC/blob/main/imgs/gradient_img.png" width="80%" />
</p>


## Fitting configurations
In our paper we fit two different simulation-informed signal models. We refer to these as two different _fitting configurations_. These correspond to forward models 1 and forward 2, namely:

- in forward model 1 (fitting configuration number 9 in our code) we estimate `fin, vCS, D0in, D0ex`;

- in forward model 2 (fitting configuration number 14 in our code) we estimate `fin, vCS, D0in, D0ex, kappa`.

Above, the parameters have the following meaning:
* f<sub>in</sub>: intra-cellular signal fraction
* vCS: characteristic volume-weighted cell diameter
* D<sub>0|in</sub>: intrinsic intra-cellular diffusivity
* D<sub>0|ex</sub>: intrinsic extra-cellular diffusivity.
* $\kappa$: cell membrane permeability

The parameters in forward model 1 were chosen so that we could compare the performance of simulation-informed fitting vs a well-established two-compartment model, accounting for restricted diffusion within cylinders and hindered extra-cellular Gaussian diffusion. Remember that our substrates where obtained on 2D histology, and thus feature cylndrical symmetry (that is why we used cylinders, instead of spheres). Note that the analytical model can only be fitted on PGSE data, and not on TRSE.

Forward model 2 instead investigates the feasibility of data-driven estimation of the microstructural parameters while also estimating the cell membrane permeability. As an example here we show how to fit `forward model 2` using the `PGSEin` protocol.


### Monte Carlo simulation-informed fitting
We provide the code to perform model fitting on the simulated signals. 

The script [`run_dictionary_fit.py`](https://github.com/radiomicsgroup/dMRIMC/blob/main/parameter_estimation/run_dictionary_fit.py) performs the MC-informed fitting for all cases. This script relies on the [`mri2micro_dictml.py`](https://github.com/radiomicsgroup/dMRIMC/blob/main/parameter_estimation/mri2micro_dictml.py) tool, released as part of [bodymritools](https://github.com/fragrussu/bodymritools) (script [mri2micro_dictml.py](https://github.com/fragrussu/bodymritools/blob/main/mrifittools/mri2micro_dictml.py)). Note that _mri2micro_dictml.py_ can be used to fit **any equation-free, numerical signal model, given examples of signals and corresponding tissue parameters for a given acquisition protocol**. 

To run `run_dictionary_fit.py`, simply clone this repository, navigate to [`parameter_estimation/`](https://github.com/radiomicsgroup/dMRIMC/tree/main/parameter_estimation/run_dictionary_fit.py) and run it:

```
python run_dictionary_fit.py
```

### Plotting fitting results
We also provide a script to generate scatter density plots that correlate estimated vs ground truth tissue parameters:

- [`density_plots_and_correlations.py`](https://github.com/radiomicsgroup/dMRIMC/blob/main/parameter_estimation/density_plots_and_correlations.py): results from MC-informed fitting for a desired protocol, SNR and forward model. You can use it like this:
```
python density_plots_and_correlations.py

```

With this, you can generate Fig.4, panels (a)-(e) related to the _in silico_ parameter estimation experiments of our paper.

<div align="center">
    <img src="https://github.com/radiomicsgroup/dMRIMC/blob/main/parameter_estimation/figs/PGSEin_SNR50_config14_all_params.jpg" alt="mcinformed" width="auto" height="auto">
</div>


## Dependencies
The code was developed with the following package versions:
- **python** (3.10.13)
- **numpy** (1.26)
- **scipy** (1.11.3)
- **seaborn** (0.12)
- **nibabel** (5.1.0)
- **mrtrix3** (3.0.4)
- **numba** (0.60.0)