## Informed microstructure parameter estimation
Signals generated from the manually segmented substrates can be used to inform parameter estimation. This is a detailed tutorial on how to use the uploaded code, as well as a way to replicate some of our figures that appear in the paper.

### Data
Since we want also to reproduce the figures of the paper we will be using real data from patients which we cannot share, so the tutorial will assume that we have generated signals from our substrates, and created niftii files partitioning the data according to the leave-one-out method, both for the generated signals and the corresponding parameters. The signals are also assumed to be corrupted with Rician noise (citation) and noise maps have been generated for each niftii using `dwidenoise` from [MRtrix3](https://mrtrix.readthedocs.io/en/latest/reference/commands/dwidenoise.html).

We have two parameter configurations corresponding to forward model 1 and forward model 2 in our paper:

Forward model 1 parameter configuration: $`f_{in},`$ vCS, $`D_{0|in}, D_{0|ex}`$

Forward model 2 parameter configuration: `$f_in$, mCS, varCS, skewCS, $D_{0|in}$, $D_{0|ex}$` 
