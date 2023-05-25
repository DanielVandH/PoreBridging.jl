# PoreBridging

This repository contains the code and data for VandenHeuvel et al. (2023). The abstract for this paper is:

> Understanding how tissue growth in porous scaffolds is influenced by geometry is a fundamental challenge in the field of tissue engineering. We investigate the influence of pore geometry on tissue growth using osteoblastic cells in 3D printed melt electrowritten scaffolds with square-shaped pores and non-square pores with wave-shaped boundaries. Using a reaction-diffusion model together with a likelihood-based uncertainty quantification framework, we quantify how the mechanisms of cell migration and cell proliferation drive tissue growth for each pore geometry. Our results show that the rates of cell migration and cell proliferation appear to be independent of the pore geometries considered, suggesting that simple square-shaped pores can be used to estimate parameters and make predictions about tissue growth in more realistic pores with complicated shapes. Our findings have important implications for the development of predictive tools for tissue engineering and experimental design, highlighting new avenues for future research.

The following sections describe:

1. The repository structure.
2. Steps for reproducing the results.

Moreover, this repository makes heavy use of my other three packages:

- [DelaunayTriangulation.jl](https://github.com/DanielVandH/DelaunayTriangulation.jl): Mesh generation.
- [FiniteVolumeMethod.jl](https://github.com/DanielVandH/FiniteVolumeMethod.jl): Solving the PDEs.
- [ProfileLikelihood.jl](https://github.com/DanielVandH/ProfileLikelihood.jl): Performing the profile likelihood analysis.

If there are any problems with the instructions that follow, please file an issue. (Issues relating to the use of the packages mentioned above should be given in those repositories, not here.)

# Repository Structure 

The repository is broken into three folders:

- `data/`

The data folder stores four files:

-- `dataset.csv`: This is an Excel Spreadsheet that stores the summary statistics for each pore considered. The `fileName` column stores references to files in folders that unfortunately cannot be shared on this repository. This spreadsheet is not so used throughout the repository, as we recompute the summary statistics directly from the identified boundaries from the MATLAB script `process_images.m` anyway, but it is useful for splitting up the data by geometry.
-- `process_images.m`: This is the MATLAB script that processes all the images. You will not be able to run this directly as the images are unfortunately not on this repository, but it serves as a useful reference. As with `dataset.csv`, the summary statistics computed in this script are not used, but the boundary representations are. Unfortunately, these boundary representations are not provided here also.
-- `square_pore_data.jld2`: This is a `.jld2` file, which you can read using [JLD2.jl](https://github.com/JuliaIO/JLD2.jl). It stores the summary statistics for the square pores considered in the paper.
-- `wave_pore_data.jld2`: Similar to `square_pore_data.jld2`, except for the wave geometry.

The `.jld2` files are the results of applying the `get_pore_data` function (available in the `src` folder) to the original raw data.

- `src/`

This folder contains the actual module `PoreBridgingProfiling` that we use for performing the analyses in the paper. Not all functions are used due to some changes we made throughout writing the paper, so the best way to interpret what we do would be to run through the actual analysis in the `paper` folder (described next), or see the exported functions in the module file. Regardless, this folder stores multiple files:

-- `PoreBridgingProfiling.jl`: This is the definition of the module that includes all the dependencies and includes all the necessary files.
-- `structs.jl`: This defines the structs that we use for storing the data and running the analyses.
-- `profile_analysis_functions.jl`: Some functions for performing the profile analysis. There is some overlap here with the full likelihood analysis.
-- `full_likelihood_analysis_functions.jl`: Some functions for performing the full likelihood analysis.
-- `data_setup.jl`: Sets up the data to be used in the likelihood function.
-- `geometry_setup.jl`: Sets up the geometry to be used for either the square, wave, or cross pore geometries.
-- `mle.jl`: Functions for computing the MLEs.
-- `objective_function.jl`: Functions for computing the log-likelihood.
-- `optimisation_setup.jl`: Functions for defining the likelihood problem.
-- `pde_setup.jl`: Functions for defining the PDE problem.
-- `prediction_intervals.jl`: Functions for computing the prediction intervals of the summary statistics, and the probability density functions for the bridging time.
-- `profile.jl`: Functions for computing the profile log-likelihoods.
-- `snapshots.jl`: Functions for computing the tissue growth variability from the profile likelihood results.
-- `summary_statistics.jl`: Functions for computing summary statistics from the numerical simulations.

- `paper/`

This is the folder that we use for actually obtaining the results in the paper. The main file to consider is

-- `complete_script.jl`: This script runs all the other scripts. It has a commented line `include(joinpath(script_path, "profile_results.jl"))`, disabling the automatic running of the profile likelihood analysis script that takes an _EXTREMELY_ long time to run. If you really want to run it, either uncomment this line and rerun `complete_script.jl` or go into the file itself and run the code. After this line, four plotting scripts are run, using [JLD2.jl](https://github.com/JuliaIO/JLD2.jl) to run the saved results from `profile_results.jl` provided to you in the `paper/jld2_files/` folder. Last, the hypothetical study performed in the paper is run from `hypothetical_study.jl`.

Now we describe all the other scripts.

-- `profile_analysis.jl`: As mentioned, this script performs the profile likelihood analyses and saves all the results as `.jld2` files using [JLD2.jl](https://github.com/JuliaIO/JLD2.jl). The script takes an extremely long time to complete. It starts by obtaining all the results with the full likelihood, and then obtains the profile likelihood results.
-- `plotting_schematics.jl`: This script gets the figure in the paper that shows the schematics for the geometries and meshes used.
-- `plotting_likelihood_analysis_results.jl`: By loading in the `.jld2` files saved from `profile_analysis.jl`, we plot in this script all the results for the likelihood analysis.
-- `plotting_reparametrisation_discussion.jl`: This script shows the difference between the $(D, \lambda)$ and $(D\lambda,\lambda)$ parametrisations of the likelihood, showing the log-likelihood surfaces in each case.
-- `plotting_tissue_growth_predictions.jl`: This script plots the predictions for the tissue growth variability.
-- `hypothetical_study.jl`: This script is used for running the hypothetical study on the cross pore geometry given in the paper, using results from the profile likelihood analysis on the square pore. The script is hopefully sufficiently commented for you to understand, especially in case you want to extend the ideas here to another type of geometry, but if you have any problems here please either email or post an issue in the issues tab.

In addition to these scripts, we also have two folders:

-- `paper/figures/`: This doesn't actually have any content, but I've kept it here in case you run `complete_script.jl` which assumes this directory exists. Most of the figures for the paper will be saved into here after you run `complete_script.jl`. 
-- `paper/jld2_files/`: This folder stores all the `.jld2` files from the profile likelihood analysis and from the hypothetical study. Everything is hopefully self-explanatory, except in `paper/jld2_files/full_likelihood_results/` where the fle names ending in `_i` refer to $u_0 = 0.2$, $u_0 = 0.3$, and $u_0=0.4$ for `i = 1`, `i = 2`, and `i = 3`, respectively.

# Reproducing the Results 

As mentioned above, the `paper/` together with the `paper/complete_script.jl` script is used for reproducing our results. To actually run this, though, requires some extra steps. 

1. You should have Julia installed. You can download Julia [here](https://julialang.org/downloads/), preferably `v1.9.0` (the latest version at the time of writing). My preferred installer for Julia is [juliaup](https://github.com/JuliaLang/juliaup).
2. You need an actual editor for Julia, e.g. [VS Code](https://code.visualstudio.com/) with the [Julia extension](https://code.visualstudio.com/docs/languages/julia).
3. Now fork or clone this repository, and open the folder in your Julia editor. Once you have done this, activate the `PoreBridgingProfiling` environment by entering `] activate .` into the Julia REPL. If done correctly, typing `]` into the REPL (without executing anything) should show `(PoreBridgingProfiling) pkg>`.
4. With the package now activated, type `] instantiate` to install all the required dependencies for the package. You can then enter `using PoreBridgingProfiling` to test if you've instantiated everything correctly - you might need to do e.g. `] precompile`; it'll tell you what to do if you've made any errors.
5. Once this is all done, you should be able to run `complete_script.jl` as needed.

Note that this should all be perfectly reproducible thanks to the `Manifest.toml` and `Project.toml` files in the repository.
