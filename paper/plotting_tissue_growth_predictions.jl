################################################################################################
##
## SETUP PACKAGES AND PARAMETERS
##
################################################################################################
using PoreBridgingProfiling
using Setfield
using JLD2
using CairoMakie
using ProfileLikelihood
using DelaunayTriangulation

jld2_path_profile = normpath(@__DIR__, "..", "paper", "jld2_files", "profile_results")
jld2_path_figures = normpath(@__DIR__, "..", "paper", "figures", "tissue_growth")

u0 = [0.2, 0.3, 0.4]
pore_data, save_options = get_data_structs();

############################################
## Square full likelihood analysis: Area and perimeter
#############################################
square_area_perimeter_profile_results = load(joinpath(jld2_path_profile, "square_area_perimeter_results.jld2"))["square_area_perimeter_results"][1];
fig_square, fig_wave = plot_tissue_growth_predictions(square_area_perimeter_profile_results, pore_data, save_options)
save(joinpath(jld2_path_figures, "square_perimeter_area_square_tissue_growth_results.pdf"), fig_square)
save(joinpath(jld2_path_figures, "square_perimeter_area_wave_tissue_growth_results.pdf"), fig_wave)

############################################
## Square full likelihood analysis: Area only
#############################################
square_area_profile_results = load(joinpath(jld2_path_profile, "square_area_results.jld2"))["square_area_results"][1];
fig_square, fig_wave = plot_tissue_growth_predictions(square_area_profile_results, pore_data, save_options)
save(joinpath(jld2_path_figures, "square_area_square_tissue_growth_results.pdf"), fig_square)
save(joinpath(jld2_path_figures, "square_area_wave_tissue_growth_results.pdf"), fig_wave)

############################################
## Wave full likelihood analysis: Area and perimeter
#############################################
wave_area_perimeter_profile_results = load(joinpath(jld2_path_profile, "wave_area_perimeter_results.jld2"))["wave_area_perimeter_results"][1];
fig_square, fig_wave = plot_tissue_growth_predictions(wave_area_perimeter_profile_results, pore_data, save_options)
save(joinpath(jld2_path_figures, "wave_perimeter_area_square_tissue_growth_results.pdf"), fig_square)
save(joinpath(jld2_path_figures, "wave_perimeter_area_wave_tissue_growth_results.pdf"), fig_wave)

############################################
## Wave full likelihood analysis: Area only
#############################################
wave_area_profile_results = load(joinpath(jld2_path_profile, "wave_area_results.jld2"))["wave_area_results"][1];
fig_square, fig_wave = plot_tissue_growth_predictions(wave_area_profile_results, pore_data, save_options)
save(joinpath(jld2_path_figures, "wave_area_square_tissue_growth_results.pdf"), fig_square)
save(joinpath(jld2_path_figures, "wave_area_wave_tissue_growth_results.pdf"), fig_wave)

