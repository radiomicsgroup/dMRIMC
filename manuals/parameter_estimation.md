## Simulation-informed microstructure parameter estimation
Signals generated from the manually segmented substrates can be used to inform parameter estimation. This is a detailed tutorial on how to use the uploaded code, as well as a way to replicate some of our figures that appear in the paper.

**Under development! Everything will be here in the coming days -- apologies for the inconvenience!**

### Data
Since we want also to reproduce the figures of the paper we will be using real data from patients which we cannot share, so the tutorial will assume that we have generated signals from our substrates, and created niftii files partitioning the data according to the leave-one-out method, both for the generated signals and the corresponding parameters. The signals are also assumed to be corrupted with Rician noise (citation) and noise maps have been generated for each niftii using `dwidenoise` from [MRtrix3](https://mrtrix.readthedocs.io/en/latest/reference/commands/dwidenoise.html). The resulting niftii files are in `leave_one_out/niftiis` (MC generated signals) and `leave_one_out/analytical` (signals from analytical expression) and the leave-one-out-partitioned signals are stored in `leave_one_out/signal_np_files`. The corresponding parameters for each configuration (see below) are in `leave_one_out/param_np_files`. Included in the data are signals generated for all three protocols mentioned in the paper, namely `PGSEin` (in silico/in vivo), `PGSEex` (ex vivo) and `TRSE` (in silico/in vivo).


### Fitting with forward model 1 (comparison with analytical expression)
We have two parameter configurations corresponding to forward model 1 and forward model 2 in our paper:

Forward model 1 parameter configuration (config number 9): $`f_{in},`$ vCS, $`D_{0|in}, D_{0|ex}`$

Forward model 2 parameter configuration (config number 13): $`f_{in},`$ mCS, varCS, skewCS, $`D_{0|in}, D_{0|ex}`$

