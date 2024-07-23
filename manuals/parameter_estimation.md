## Simulation-informed microstructure parameter estimation
Signals generated from the manually segmented substrates can be used to inform parameter estimation. This is a detailed tutorial on how to use the uploaded code, as well as a way to replicate some of our figures that appear in the paper.

**Under development! Everything will be here in the coming days -- apologies for the inconvenience!**

### Data
Since we want also to reproduce the figures of the paper we will be using real data from patients which we cannot share, so the tutorial will assume that we have generated signals from our substrates, and created niftii files partitioning the data according to the leave-one-out method, both for the generated signals and the corresponding parameters. The signals are also assumed to be corrupted with Rician noise (citation) and noise maps have been generated for each niftii using `dwidenoise` from [MRtrix3](https://mrtrix.readthedocs.io/en/latest/reference/commands/dwidenoise.html). The resulting niftii files are in `leave_one_out/niftiis` (MC generated signals) and `leave_one_out/analytical` (signals from analytical expression) and the leave-one-out-partitioned signals are stored in `leave_one_out/signal_np_files`. The corresponding parameters for each configuration (see below) are in `leave_one_out/param_np_files`. Included in the data are signals generated for all three protocols mentioned in the paper, namely `PGSEin` (in silico/in vivo), `PGSEex` (ex vivo) and `TRSE` (in silico/in vivo).


### Fitting configurations
We have two parameter configurations corresponding to forward model 1 and forward model 2 in our paper:

- Forward model 1 parameter configuration (config number 9): $`f_{in},`$ vCS, $`D_{0|in}, D_{0|ex}`$

- Forward model 2 parameter configuration (config number 13): $`f_{in},`$ mCS, varCS, skewCS, $`D_{0|in}, D_{0|ex}`$

The parameters in forward model 1 were chosen so that we could compare the performance of numerically-obtained signals vs the analytical expression in the context of model fitting. The ones in forward model 2 were chosen to investigate the first three moments of the cell size distibution.

### Fitting procedures
1. Fitting the analytical model for $$$$$$. This script performs the fitting of the analytical expression using all the signals at the same time.

  - `python dri2mc_maxlikcyl.py --noise analytical/noise_maps/PGSEin_noise_analytical_SNR50_all_signals.nii --modstr DinDex --pmin 8.0,0.8,0.0,0.5 --pmax 20.0,3.0,0.9,3.0 --sldim 0 --nw 12 --ncpu 10 analytical/niftiis/PGSEin_all_signals_config_9_SNR_50.nii PGSEin_for_maxlik.bval res_maxlik/old_way_cylinders/config_9/SNR_50`

2. Fitting using the MC-generated signals with the PGSEin protocol. The `run_all.py` script will run the fitting routine for all protocols and configurations.

  - `python run_all.py`

To reproduce Fig 5 from our paper i.e contour plots comparing the MC-informed parameter estimation with the analyical expression for a hindered-restricted diffusion model use `explore_niftii_LeaveOneOut.py` and `explore_niftii_LeaveOneOut_ANALYTICAL.py` after the fittings.

1. `python explore_niftii_LeaveOneOut.py PGSEin 50 9`: This will create Fig 5a (MC-informed fitting)
2. `python explore_niftii_LeaveOneOut_ANALYTICAL.py PGSEin 50 9`: This will create Fig 5b (analytical fitting)
