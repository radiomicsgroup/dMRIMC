import numpy as np
import argparse

# Load combined parameter array
param_arr = np.load("all_parameters_all_substrates.npy")

# Parameters are stored in this order:
# fin, mCS, varCS, skewCS, vCS_sph, vCS_cyl, Din, Dex, kappa
param_index = {
    "fin":0,
    "mCS":1,
    "varCS":2,
    "skewCS":3,
    "vCS_sph":4,
    "vCS_cyl":5,
    "Din":6,
    "Dex":7,
    "kappa":8,
    }

# Argparse setup
parser = argparse.ArgumentParser()
parser.add_argument('--params', type=str, required=True, 
                    choices=['fin', 'mCS', 'varCS', 'skewCS', 'vCS_sph', 
                             'vCS_cyl', 'Din', 'Dex', 'kappa'],
                    nargs='+',
                    help='Parameter names to keep')
args = parser.parse_args()

p_columns = args.params
p_columns.append('Din')
p_columns.append('Dex')
p_columns.append('kappa')
cols_to_keep = [param_index[p] for p in p_columns]

param_arr_filtered = param_arr[:, cols_to_keep]

np.save('param_arr_subset.npy', param_arr_filtered)