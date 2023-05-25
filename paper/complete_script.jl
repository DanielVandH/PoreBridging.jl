script_path = normpath(@__DIR__, "..", "paper")

## Only run if you want to rerun all the profile and likelihood analyses. 
## This will take an EXTREMELY LONG TIME. 
# include(joinpath(script_path, "profile_results.jl"))

## Plotting scripts 
file_names = [
    "plotting_schematics.jl",
    "plotting_likelihood_analysis_results.jl",
    "plotting_reparametrisation_discussion.jl",
    "plotting_tissue_growth_predictions.jl",
]
for file_name in file_names
    include(joinpath(script_path, file_name))
end

## The hypothetical study 
include(joinpath(script_path, "hypothetical_study.jl"))