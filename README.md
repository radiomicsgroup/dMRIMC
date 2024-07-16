
# A Monte Carlo simulation framework for histology-informed diffusion MRI cancer characterisation and microstructural parameter estimation

<div align="center">
  <img src="https://github.com/radiomicsgroup/dMRIMC/blob/main/imgs/diagram22.png" alt="qrcode" width="auto" height="auto">
</div>


<div align="center">
  <img src="https://github.com/radiomicsgroup/dMRIMC/blob/main/imgs/qr_img_MC_2024_paper.png" alt="qrcode" width="auto" height="auto">
</div>

## [MedArxiv preprint](https://www.medrxiv.org/content/10.1101/2024.07.15.24310280v1)

# UNDER CONSTRUCTION

This code accompanies our [publication](https://www.medrxiv.org/content/10.1101/2024.07.15.24310280v1), sharing the process of creating realistic cancer substrates, running Monte Carlo simulations in them and synthesizing the resulting MRI signal according to any user-specified PGSE protocol. The process can be broken down in the following stages:
- Select region of interest (ROI) from histology (e.g using [QuPath](https://qupath.github.io/)) and take a high resolution screenshot making sure the scale bar is in view.
- Segment the features of interest using [Inkscape](https://inkscape.org/) and save the file as `svg`
- Import the `svg` into [Blender](https://www.blender.org/) and obtain geometry files: 1 file for the whole substrate for extracellular simulations and n files corresponding to each individual cell.
- Use the `ply` files with the [MCDC simulator](https://github.com/jonhrafe/MCDC_Simulator_public) to get the random walks for all objects.
- Synthesize the total signal (extracellular + intracellular) according to any PGSE protocol


The generated signals can then be used to perform inference and obtain parametric maps of scans. Note that in this repository as part of an example the segmented cells of a mouse histology sample are included in a zip file inside the `playgrounds` folder. All other files (such as trajectory files or signal files) are not included for practical reasons as the size would be too big to host here.



## License
This repository is distributed under the **Attribution-NonCommercial-ShareAlike 4.0 International license** ([CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0)). Copyright (c) 2024, Fundació Privada Institut d’Investigació Oncològica de Vall d’Hebron (Vall d'Hebron Institute of Oncology (VHIO), Barcelona, Spain). All rights reserved. Link to license [here](). 

## From histological image to MRI signal
### Segmenting with Inkscape
Inkscape is a very versatile software with its primary use being vector graphics. It is useful in our case because it can output `svg` files that store the information about their geometry in paths and can thus be converted to meshes. To perform the segmentation we used the pen tool as it allowed for the fast creation of arbitrary shapes.

Once all desired features are segmented, the image has to be cropped to contain the segmented area only. This can be done by a Clip operation, just draw a rectangle over the segmented area (making sure it is higher up in the layer hierarchy) and with both it and the background image selected go to Object-->Clip-->Set clip. To also resize the canvas of the image, select everything with the mouse and then press Ctrl+D to open the document properties menu and click "Resize to content". Save the cropped segmented image as a separate file.

Now the features of the image can be exported. We need two kinds of exports, one is for the whole features together for the extracellular simulation and a batch export of all the cells. To get the `svg` for the whole substrate open the cropped image and combine all the paths (i.e cellular features) into one path by selecting all the paths and pressing Ctrl+K (or go to Path-->Combine). With all the features now in one path open the export panel and export the path only (with the background image) as an Inkscape `svg` file with `_ALL_STRUCTURES` appended to the end. For the individual cells, open the original image (or the cropped image before combining the paths) and in the export panel select Batch Export. To be safe, select all the cells and tick the Export Selected Only box. Create a folder inside the working directory called `cells` then export all the cells as Inkscape `svg` files in a folder called `svgs`. If you have other water-bearing features (such as lumen spaces) that you want to run a simulation in, make a new folder (e.g `lumen`) inside the working directory and repeat the process.

Finally, we need the areas of the cells and the entire substrate to weigh their signals appropriately later. This is done with `make_cell_vol_array.py` and `make_extra_vol_array.py`

## From `svg` to `ply` in Blender
To create the `ply` files we use the scripts in the `geometry_scripts` folder. Before we run any scripts first we must establish the correct scaling factor for each substrate. When opening the files in Blender they will be slightly different from Inkscape so we need to apply a correction. To do that measure the same object in both Blender and Inkscape and using the scale bar from the histology image in Inkscape determine the correction factor to multiply the geometry with. With that in hand, open Blender and navigate to the scripting tab. There are two pairs of files, one for batch processing and one for single file processing (i.e for the entire substrate, the extracellular part). Both pairs need to be located in the same directory as the folder containing the substrate images and exported `svgs`. Inside each file the name of the folder containing the substrate needs to be specified, as well as the scale factor. First run the `svg_to_stl.py` scripts and then the `stl_to_ply.py` scripts. If all went well, we should have all we need for the simulations!

## Run simulations with MCDC
For each substrate you need 2 folders, one for extracellular simulations and one for the intracellular. To run simulations first create a folder with 2 subfolders: e.g  for `substrate_name_CELLS` create `substrate_name_CELLS/random_walks` and `substrate_name_CELLS/conf_files`.  The first one will house the trajectories of the spins and the latter the simulation configuration files. For more information on MCDC check the repository (link). For the cells (and any other intracellular geometry) run `create_config_files_MISC.py` (see example). For the extracellular run `create_config_files_ALL_STRUCTURES.py`. Inside both scripts you can specify the parameters needed by MCDC (diffusivities, simulation duration etc).

With the configuration files ready the simulations can be run using the `run_sims_in_parallel.py` script.  It requires the destination folder's name, the substrate type that we want (`intra` or `extra`) and the number of cores that can be used. The script will try to run one simulation per core using a queue if the simulations asked for are more than the number of cores. It will print out information while running. When the simulations are done you will see the number of errors (if any) and then you can run the script `manual_testing.py` to scan the resulting files for inconsistencies (check the file for info on the tests).

We now have the trajectories of the spins in the `.traj` files inside the `random_walks` folders.

Code to recreate the mouse example contained in the repo:

Create configuration files for the simulations

`python create_config_files_MISC.py mouse_example_CELLS/ mouse_example cells`

`python create_config_files_ALL_STRUCTURES.py mouse_example_EXTRA/ mouse_example`


Run simulations (both intra- and extracellular) using 10 threads

`python run_sims_in_parallel.py mouse_example_CELLS intra 10`

`python run_sims_in_parallel.py mouse_example_EXTRA extra 10`


Run some sanity checks (optional)

`python manual_testing.py mouse_example_EXTRA intra`

`python manual_testing.py mouse_example_EXTRA extra`

## Signal synthesis
The scripts included can synthesize any PGSE signal provided the parameters. Create a folder `sequence parameters` and in there a subfolder for any PGSE protocol you want to simulate. The subfolder has to contain four text files: 
- `custom.bval` - bvalues ($s/mm^2$)
- `custom.gdur1` - gradient duration 1 ($ms$)
- `custom.gdur2` - gradient duration 2 ($ms$)
- `custom.gsep12` - gradient separation 1-2 ($ms$)
Each file should contain the parameter values for each measurement separated with a space. To perform the synthesis, run `run_all_for_misc.py` where misc stands for misc protocol. It requires the name of the protocol corresponding to the name of the folder where the files describing said protocol are (bvalue, timings), plus the substrate type (`cells` or `EXTRA`). It assumes that the trajectories' folder is located in `SIMULATIONS`. By default `run_all_for_misc.py` will look at the `playgrounds` folder and synthesize signals for all simulations matching the substrates contained there.

After synthesizing both the intracellular and extracellular signal, use `aggregate_all_intra.py` and `aggregate_all_extra.py` to volume weight all the signals and collect them in one file for intracellular and one for extracellular. Note that this requires having the area of each cellular feature (e.g cells) in a dictionary for reading during runtime (for an example check `vol_INTRA_func_CUSTOM.py`).

For the mouse example:

Synthesize signals for the `mouse_example` substrate according to the `CUSTOM_PGSE` protocol for `EXTRA`cellular diffusion

`python run_all_for_misc.py CUSTOM_PGSE EXTRA`

Same but for intracellular diffusion (cells)

`python run_all_for_misc.py CUSTOM_PGSE cells`

Aggregate all extracellular and intracellular signals

`python aggregate_all_extra.py CUSTOM_PGSE CUSTOM_PGSE`

`python aggregate_all_intra.py CUSTOM_PGSE CUSTOM_PGSE`

Now the in `SIGNAL_SYNTHESIS/aggregated_signals` we have two folders containing the intra- and extracellular components of the signal for each diffusivity value. Combining them with appropriate area/volume-weighting gives us the total measured signal of the simulated voxel.

## Code
### Simulations
- `create_config_files_ALL_STRUCTURES.py`: Creates MCDC configuration files for the whole substrate (extracellular simulation)
- `create_config_files_MISC.py`: Creates MCDC configuration files for each cell in the substrate, or for each element of a cellular feature such as lumen (intracellular simulation)
- `manual_testing.py`: Runs some post-simulation tests on the trajectory files
- `parallel_execution_functions.py`: Contains three variations of the function to run the simulations in parallel: Just running them in parallel, running them with a time limit in case they get stuck and retrying those that got stuck, running them with a time limit but with no automatic retry
- `post_simulation_tests.py`: Contains the post simulation testing code
- `retry_all_sick.py`: Goes through the error logs of a simulation and reruns the ones that got stuck the first time
- `run_sims_in_parallel.py`: Main simulation running script, read the `conf_files` folder and calls the parallel execution function.
### Signal synthesis
- `aggregate_all_extra.py`:  Aggregates all extracellular features, can be only extracellular space simulation and should also be modified to include other features that are water-bearing but are considered part of the extracellular signal (e.g lumen)
- `aggregate_all_intra.py`: Aggregates all intracellular features
- `MRI_functions.py`: Contains the functions:
	- `MRI_functions.load_trajectories_2D`: Loads the trajectories from MCDC and prepares them in a convenient way
	- `MRI_functions.construct_signal_PGSE`: Simulates the MR signal from the set of random walks coming from the trajectory files and returns MRI signal at infinite SNR
	- `MRI_functions.synthesis_misc_PGSE`: Function calling `MRI_functions.construct_signal_PGSE`
- `run_all_for_misc.py`: Runs the synthesis process, asking for the protocol name as it appears in the folder name with the parameters (bval etc) and the substrate type `cells` or `EXTRA`. Will also keep track of any errors.
- `synthesis_PGSE_one_file.py`: Standalone function to 
- `vol_EXTRA_func_CUSTOM.py`: Contains the code for extracellular aggregation
- `vol_INTRA_func_CUSTOM.py`: Contains the code for intracellular aggregation
### Geometric operations
- `svg_to_2D_stl.py`:  Imports an `svg`, converts it to an `stl` and exports it
- `stl_to_ply.py`: Imports `stls` and converts them to `ply` files. Going directly from `svg` to `ply` is not robust
- `svg_to_2D_stl_SIGNLE_FILE.py`: Same as above but for a single file, used for the big extracellular space
- `stl_to_ply_SIGNLE_FILE.py`: Same as above but for a single file, used for the big extracellular space
- `make_cell_vol_array.py`: Goes through the cell `ply` files and creates a dictionary with their areas to be used for signal weighting during aggregation (see above)
- `make_extra_vol_array.py`: Same as above but for the extracellular space

## Dependencies
The code was developed with the following package versions:
- **python** (3.9.17)
- **numpy** (1.25.2)
- **scipy** (1.10.1)
- **pyntcloud** (0.3.1)
- **opencv** (4.7.0)

Blender version `3.6.12` was used for the geometric operations

