## From histological images to diffusion MRI signals
The process can be broken down in the following stages:

- Select region of interest (ROI) from histology (e.g using [QuPath](https://qupath.github.io/)) and take a high resolution screenshot making sure the scale bar is in view.
- Segment the features of interest using [Inkscape](https://inkscape.org/) and save the file as `svg`
- Import the `svg` into [Blender](https://www.blender.org/) and obtain geometry files: 1 file for the whole substrate for extracellular simulations and n files corresponding to each individual cell.
- Use the `ply` files with the [MCDC simulator](https://github.com/jonhrafe/MCDC_Simulator_public) to get the random walks for all objects.
- Synthesize the total signal (extracellular + intracellular) according to any PGSE protocol


The generated signals can then be used to perform inference and obtain parametric maps of scans. Note that in this repository as part of an example the segmented cells of a mouse histology sample are included inside the `playgrounds` folder. All other files (such as trajectory files or signal files) are not included for practical reasons as the size would be too big to host here.

### Segmenting with Inkscape
<p align="center">
  <img src="https://github.com/radiomicsgroup/dMRIMC/blob/main/imgs/ex.png" width="45%" />
  <img src="https://github.com/radiomicsgroup/dMRIMC/blob/main/imgs/ex2.png" width="45%" />
</p>

Inkscape is a very versatile software with its primary use being vector graphics. It is useful in our case because it can output `svg` files that store the information about their geometry in paths and can thus be converted to meshes. To perform the segmentation we used the pen tool as it allowed for the fast creation of arbitrary shapes.

Once all desired features are segmented, the image has to be cropped to contain the segmented area only. This can be done by a Clip operation, just draw a rectangle over the segmented area (making sure it is higher up in the layer hierarchy) and with both it and the background image selected go to Object-->Clip-->Set clip. To also resize the canvas of the image, select everything with the mouse and then press Ctrl+D to open the document properties menu and click "Resize to content". Save the cropped segmented image as a separate file.

Now the features of the image can be exported. We need two kinds of exports, one is for the whole features together for the extracellular simulation and a batch export of all the cells. To get the `svg` for the whole substrate open the cropped image and combine all the paths (i.e cellular features) into one path by selecting all the paths and pressing Ctrl+K (or go to Path-->Combine). With all the features now in one path open the export panel and export the path only (with the background image) as an Inkscape `svg` file with `_ALL_STRUCTURES` appended to the end. For the individual cells, open the original image (or the cropped image before combining the paths) and in the export panel select Batch Export. To be safe, select all the cells and tick the Export Selected Only box. Create a folder inside the working directory called `cells` then export all the cells as Inkscape `svg` files in a folder called `svgs`.

Exporting all elements individually is necessary to calculate the metrics of each substrate (e.g average cell size, intracellular fraction etc). Running [make_cell_vol_array.py](https://github.com/radiomicsgroup/dMRIMC/blob/main/geometry_scripts/make_cell_vol_array.py) and [make_extra_vol_array.py](https://github.com/radiomicsgroup/dMRIMC/blob/main/geometry_scripts/make_extra_vol_array.py) gives a file with the areas of all cells and the area of the extracellular space respectively. Then we can calculate various metrics for the substrate by running `get_metrics.py` (TBA). This file is necessary to establish the ground truth for the parameter estimation and create the parameters' dictionary. In this repository a parameters' dictionary is already provided created from the substrates discussed in the paper. The substrates themselves have also been released on [Zenodo](https://doi.org/10.5281/zenodo.14559103)

## From `svg` to `ply` in Blender
To create the `ply` files we use the scripts in the `geometry_scripts` folder. Before we run any scripts first we must establish the correct scaling factor for each substrate. When opening the files in Blender they will be slightly different from Inkscape so we need to apply a correction. To do that measure the same object in both Blender and Inkscape and using the scale bar from the histology image in Inkscape determine the correction factor to multiply the geometry with. With that in hand, open Blender and navigate to the scripting tab. There are two pairs of files, one for batch processing and one for single file processing (i.e for the entire substrate, the extracellular part). Both pairs need to be located in the same directory as the folder containing the substrate images and exported `svgs`. Inside each file the name of the folder containing the substrate needs to be specified, as well as the scale factor. First run the [svg_to_stl.py](https://github.com/radiomicsgroup/dMRIMC/blob/main/geometry_scripts/svg_to_stl.py) scripts and then the [stl_to_ply.py](https://github.com/radiomicsgroup/dMRIMC/blob/main/geometry_scripts/stl_to_ply.py) scripts. If all went well, we should have all we need for the simulations!

## Run simulations with MCDC
<p align="center">
  <img src="https://github.com/radiomicsgroup/dMRIMC/blob/main/imgs/mc_walkers.gif" width="60%" />
</p>


To run simulations we first need to make the simulation configuration files (for more info check the MCDC [documentation](https://github.com/jonhrafe/MCDC_Simulator_public)), using `create_config_files.py`. The script assumes you have a folder called `simulations` to house the simulation results and a folder with your ply files `ply_files`. It will make all necessary folders and configuration files for each substrate with varying intracellular/extracellular intrinsic diffusivities and permeability ($\kappa$) values. The default configuration is:

- $D_{0|in} / D_{0|ex} \in [0.8, 3] \textmu m^2/ms$
- $\kappa \in [0, 40] \textmu m/s$ ([more info here](https://doi.org/10.1002/mrm.29720))

With the configuration files ready the simulations can be run using the [run_sims_in_parallel.py](https://github.com/radiomicsgroup/dMRIMC/blob/main/geometry_scripts/run_sims_in_parallel.py) script. It requires the location of the MCDC binary, the destination folder's name, the substrate type that we want (`intra` or `extra`) and the number of cores that can be used. The script will try to run one simulation per core using a queue if the simulations asked for are more than the number of cores. It will print out information while running.

We now have the trajectories of the spins in the `.traj` files inside the `random_walks` folders.

Code to recreate the mouse example contained in the repo:

Create configuration files for the simulations

`python create_config_files.py`

Run the simulations using 10 threads

`python run_sims_in_parallel.py mouse_example 10`


## Signal synthesis

<div align="center">
    <img src="https://github.com/radiomicsgroup/dMRIMC/blob/main/imgs/d15.jpg" alt="signals" width="50%" height="50%">
</div>

The scripts included can synthesize any PGSE signal provided the parameters. Create a folder `sequence parameters` and in there a subfolder for any PGSE protocol you want to simulate. The subfolder has to contain four text files: 
- `custom.bval` - bvalues (s/mm<sup>2</sup>)
- `custom.gdur1` - gradient duration 1 (ms)
- `custom.gdur2` - gradient duration 2 (ms)
- `custom.gsep12` - gradient separation 1-2 (ms)

_Note: the names gdur1 and gdur2 are an artifact from previous versions of the code, they are usually the same file_

Each file should contain the parameter values for each measurement separated with a space. To perform the synthesis, run [run_synths_in_parallel.py](https://github.com/radiomicsgroup/dMRIMC/blob/main/geometry_scripts/run_all_for_misc.py). Like the script handling the simulations it can work with multiple threads to speed up processing, for more info check [standalone_synth_PGSE_for_parallel_running.py](). To run, it requires the name of the protocol corresponding to the name of the folder where the files describing said protocol are (bvalue, timings) and the number of threads to use. It assumes that the trajectories' folder is located in `simulations`. By default the script needs a list of substrates which is defined inside the scipt (here only `mouse_example`). Since we are dealing with 2D substrates, we are averaging the x,y directions

For the mouse example:

Synthesize signals for the `mouse_example` substrate according to the `CUSTOM_PGSE` protocol with 10 threads

`python run_synths_in_parallel.py CUSTOM_PGSE 10`

The synthesized signals are saved at `synthesized_signals/CUSTOM_PGSE`

With the signals synthesized and the metrics of the substrates calculated 

## Code
### Simulations
(For batch processing, it is recommended to recompile MCDC without the notification functionality)
- `create_config_files.py`: Creates MCDC configuration files for the whole substrate (extracellular simulation)
- `parallel_execution_functions.py`: Contains three variations of the function to run the simulations in parallel: Just running them in parallel, running them with a time limit in case they get stuck and retrying those that got stuck, running them with a time limit but with no automatic retry
- `run_sims_in_parallel.py`: Main simulation running script, read the `conf_files` folder and calls the parallel execution function.
### Signal synthesis
- `run_synths_in_parallel.py`: Runs the synthesis process, asking for the protocol name as it appears in the folder name with the parameters (bval etc) and the number of cores. Will also keep track of any errors.
### Geometric operations
- `svg_to_2D_stl.py`:  Imports an `svg`, converts it to an `stl` and exports it
- `stl_to_ply.py`: Imports `stls` and converts them to `ply` files. Going directly from `svg` to `ply` is not robust
- `svg_to_2D_stl_SIGNLE_FILE.py`: Same as above but for a single file, used for the big extracellular space
- `stl_to_ply_SIGNLE_FILE.py`: Same as above but for a single file, used for the big extracellular space
- `make_cell_vol_array.py`: Goes through the cell `ply` files and creates a dictionary with their areas
- `make_extra_vol_array.py`: Same as above but for the extracellular space


## Dependencies
The code was developed with the following package versions:
- **python** (3.9.17)
- **numpy** (1.25.2)
- **scipy** (1.10.1)
- **pyntcloud** (0.3.1)
- **tqdm** (4.65.0)
- **Blender** `3.6.12`
- **QuPath** `0.5`
