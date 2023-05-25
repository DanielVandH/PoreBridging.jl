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

jld2_path_full_likelihood = normpath(@__DIR__, "..", "paper", "jld2_files", "full_likelihood_results")
jld2_path_profile = normpath(@__DIR__, "..", "paper", "jld2_files", "profile_results")
jld2_path_figures = normpath(@__DIR__, "..", "paper", "figures", "likelihood_analysis")

u0 = [0.2, 0.3, 0.4]
pore_data, save_options = get_data_structs();

############################################
## Square full likelihood analysis: Area and perimeter
#############################################
square_perimeter_area_full_likelihood_results_1 = load(joinpath(jld2_path_full_likelihood, "square_perimeter_area_full_likelihood_results_1.jld2"))["square_perimeter_area_full_likelihood_results"];
square_perimeter_area_full_likelihood_results_2 = load(joinpath(jld2_path_full_likelihood, "square_perimeter_area_full_likelihood_results_2.jld2"))["square_perimeter_area_full_likelihood_results"];
square_perimeter_area_full_likelihood_results_3 = load(joinpath(jld2_path_full_likelihood, "square_perimeter_area_full_likelihood_results_3.jld2"))["square_perimeter_area_full_likelihood_results"];
all_results = [
    square_perimeter_area_full_likelihood_results_1,
    square_perimeter_area_full_likelihood_results_2,
    square_perimeter_area_full_likelihood_results_3
];
square_area_perimeter_profile_results = load(joinpath(jld2_path_profile, "square_area_perimeter_results.jld2"))["square_area_perimeter_results"];
fig = plot_square_likelihood_analysis(all_results, u0, pore_data, save_options, square_area_perimeter_profile_results)
save(joinpath(jld2_path_figures, "combined_square_area_perimeter_results.pdf"), fig)

############################################
## Square full likelihood analysis: Area only
#############################################
square_area_full_likelihood_results_1 = load(joinpath(jld2_path_full_likelihood, "square_area_full_likelihood_results_1.jld2"))["square_area_full_likelihood_results"];
square_area_full_likelihood_results_2 = load(joinpath(jld2_path_full_likelihood, "square_area_full_likelihood_results_2.jld2"))["square_area_full_likelihood_results"];
square_area_full_likelihood_results_3 = load(joinpath(jld2_path_full_likelihood, "square_area_full_likelihood_results_3.jld2"))["square_area_full_likelihood_results"];
all_results = [
    square_area_full_likelihood_results_1,
    square_area_full_likelihood_results_2,
    square_area_full_likelihood_results_3
];
square_area_profile_results = load(joinpath(jld2_path_profile, "square_area_results.jld2"))["square_area_results"];
fig = plot_square_likelihood_analysis(all_results, u0, pore_data, save_options, square_area_profile_results,
    square_ptb_xticks=(18:6:42, [L"%$s" for s in 18:6:42]),
    wave_ptb_xticks=(15:4:31, [L"%$s" for s in 15:4:31]))
save(joinpath(jld2_path_figures, "combined_square_area_results.pdf"), fig)

############################################
## Wave full likelihood analysis: Area and perimeter
#############################################
wave_perimeter_area_full_likelihood_results_1 = load(joinpath(jld2_path_full_likelihood, "wave_perimeter_area_full_likelihood_results_1.jld2"))["wave_perimeter_area_full_likelihood_results"];
wave_perimeter_area_full_likelihood_results_2 = load(joinpath(jld2_path_full_likelihood, "wave_perimeter_area_full_likelihood_results_2.jld2"))["wave_perimeter_area_full_likelihood_results"];
wave_perimeter_area_full_likelihood_results_3 = load(joinpath(jld2_path_full_likelihood, "wave_perimeter_area_full_likelihood_results_3.jld2"))["wave_perimeter_area_full_likelihood_results"];
all_results = [
    wave_perimeter_area_full_likelihood_results_1,
    wave_perimeter_area_full_likelihood_results_2,
    wave_perimeter_area_full_likelihood_results_3
];
wave_area_perimeter_profile_results = load(joinpath(jld2_path_profile, "wave_area_perimeter_results.jld2"))["wave_area_perimeter_results"];
fig = plot_wave_likelihood_analysis(all_results, u0, pore_data, save_options, wave_area_perimeter_profile_results)
save(joinpath(jld2_path_figures, "combined_wave_area_perimeter_results.pdf"), fig)

############################################
## Wave full likelihood analysis: Area only
############################################
wave_area_full_likelihood_results_1 = load(joinpath(jld2_path_full_likelihood, "wave_area_full_likelihood_results_1.jld2"))["wave_area_full_likelihood_results"];
wave_area_full_likelihood_results_2 = load(joinpath(jld2_path_full_likelihood, "wave_area_full_likelihood_results_2.jld2"))["wave_area_full_likelihood_results"];
wave_area_full_likelihood_results_3 = load(joinpath(jld2_path_full_likelihood, "wave_area_full_likelihood_results_3.jld2"))["wave_area_full_likelihood_results"];
all_results = [
    wave_area_full_likelihood_results_1,
    wave_area_full_likelihood_results_2,
    wave_area_full_likelihood_results_3
];
wave_area_profile_results = load(joinpath(jld2_path_profile, "wave_area_results.jld2"))["wave_area_results"];
fig = plot_wave_likelihood_analysis(all_results, u0, pore_data, save_options, wave_area_profile_results;
    Dλ_lim=(0.0, 450.0),
    square_ptb_xticks=(12:6:50, [L"%$s" for s in 12:6:50]),
    wave_ptb_xticks=(12:6:40, [L"%$s" for s in 12:6:40]))
save(joinpath(jld2_path_figures, "combined_wave_area_results.pdf"), fig)