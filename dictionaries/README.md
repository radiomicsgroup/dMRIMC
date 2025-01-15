## Dictionary files
This folder contains the dictionary files that were used in our [preprint](https://doi.org/10.1101/2024.07.15.24310280).

- `signal_dictionaries`: Contains .npy files that hold the signals generated for all substrates, for each protocol used in the paper: PGSEin, TRSE, PGSEex. PGSEin and TRSE have a shape of (1800, 21) and PGSEex (1152, 8).
- `parameter_dictionaries`: Contains .npy files with shapes of either (1800, 4) for forward model 1 or (1800, 6) for forward model 2. The parameters estimated are (in that order):
	- `forward_model_1`: (f<sub>in</sub>, vCS<sub>cyl</sub>, D<sub>0|in</sub>, D<sub>0|ex</sub>)
	- `forward_model_2`: (f<sub>in</sub>, mCS, vCS, skewCS, D<sub>0|in</sub>, D<sub>0|ex</sub>)

### Parameter
- `f<sub>in</sub>`: Intracellular fraction
	- `Units: n/a`
	- `Range: [0.023, 0.867]`
- `vCS_cyl`: Volume-weighted cell size for a cylindrical geometry system
	- `Units: μm`
	- `Range: [8.2, 19.9]`
- `mCS`: Mean of the cell size distribution
	- `Units: μm`
	- `Range: [6,1, 15.9]`
- `varCS`: Variance of the cell size distribution
	- `Units: μm2`
	- `Range:[2.3, 19.7]`
- `skewCS`: Skewness of the cell size distribution
	- `Units: n/a`
	- `Range:[-0.528, 0.861]`
- `D0in`
	- `Units: μm2/ms`
	- `Range: [0.8, 3.0] (PGSEin, TRSE) | [0.8, 2.5] (PGSEex)`
- `D0ex`
	- `Units: μm2/ms`
	- `Range: [0.8, 3.0] (PGSEin, TRSE) | [0.8, 2.5] (PGSEex)`

## Protocol details
- `PGSEin` - 21 measurements:
	- `bvalues (s/mm2)`: 0 50 100 400 900 1200 1500 0 50 100 400 900 1200 1500 0 50 100 400 900 1200 1500
	- `δ (ms)`: 0.000 3.900 5.200 9.200 15.000 18.200 21.000 0.000 3.900 5.200 9.200 13.000 15.800 18.500 0.000 3.900 5.200 9.200 12.900 15.800 18.400
	- `Δ (ms)`: 0.000 27.800 29.000 33.000 28.700 31.800 34.700 0.000 27.800 29.000 33.000 37.000 39.600 42.300 0.000 27.800 29.000 33.000 36.700 39.600 42.300
- `TRSE` - 21 measurements:
	- `bvalues (s/mm2)`: 0 50 100 400 900 1200 1600 0 50 100 400 900 1200 1600 0 50 100 400 900 1200 1600
	- `δ1 (ms)`: 0 8.9 8.9 8.9 8.9 8.9 8.9 0.0 13.2 13.2 13.2 13.2 13.2 13.2 0.0 18.9 18.9 18.9 18.9 18.9 18.9
	- `δ2 (ms)`: 0.0 17.6 17.6 17.6 17.6 17.6 17.6 0.0 19.3 19.3 19.3 19.3 19.3 19.3 0.0 21.0 21.0 21.0 21.0 21.0 21.0
	- `δ3 (ms)`: 0.0 20.4 20.4 20.4 20.4 20.4 20.4 0.0 24.8 24.8 24.8 24.8 24.8 24.8 0.0 30.5 30.5 30.5 30.5 30.5 30.5
	- `δ4 (ms)`: 0.0 6.0 6.0 6.0 6.0 6.0 6.0 0.0 7.7 7.7 7.7 7.7 7.7 7.7 0.0 9.5 9.5 9.5 9.5 9.5 9.5
	- `Δ12 (ms)`: 16.5 16.5 16.5 16.5 37.0 37.0 37.0 37.0
	- `Δ14 (ms)`: 16.5 16.5 16.5 16.5 37.0 37.0 37.0 37.0
- `PGSEex` - 8 measurements:
	- `bvalues (s/mm2)`: 7.94 4628.46 520.07 2063.00 7.58 4618.79 516.84 2056.56
	- `δ (ms)`: 12.0 12.0 12.0 12.0 12.0 12.0 12.0 12.0
	- `Δ (ms)`: 16.5 16.5 16.5 16.5 37.0 37.0 37.0 37.0

--

