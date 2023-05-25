################################################################################################
##
## SETUP PACKAGES
##
################################################################################################
using PoreBridgingProfiling
using Setfield
using JLD2
using CairoMakie
using ProfileLikelihood

################################################################################################
##
## FUNCTIONS FOR SAVING
##
################################################################################################
jld2_path_profiles = normpath(@__DIR__, "..", "paper", "jld2_files", "profile_results")
jld2_path_full_likelihood = normpath(@__DIR__, "..", "paper", "jld2_files", "full_likelihood_results")

function clear_results(results)
    _D = results.C_D
    _λ = results.C_λ
    _u₀ = results.C_u₀
    results = results.Cprof[:Dλ]
    θ1 = get_parameter_values(results)
    ℓ1 = get_profile_values(results)
    mle_val = results.parent.likelihood_solution[:Dλ]
    ci = get_confidence_intervals(results)
    lower_ci = first(ci)
    upper_ci = last(ci)
    return (θ=θ1, ℓ=ℓ1, mle_val=mle_val, lower_ci=lower_ci, upper_ci=upper_ci, C_D=_D, C_λ=_λ, C_u₀=_u₀)
end

################################################################################################
##
## ANALYSIS SETUP
##
################################################################################################
lb = [1e-6, 1e-6]
ub_product = [500.0, 5.0]
ub_separate = [1000.0, 5.0]
regular_grid_product = RegularGrid(lb, ub_product, 40)
regular_grid_separate = RegularGrid(lb, ub_separate, 40)
wave_ub_product = [600.0, 5.0]
wave_ub_separate = [600.0, 5.0]
wave_regular_grid_product = RegularGrid(lb, wave_ub_product, 30)
wave_regular_grid_separate = RegularGrid(lb, wave_ub_separate, 30)
u0 = [0.2, 0.3, 0.4]

pore_data, save_options = get_data_structs();
optimisation_options = OptimisationOptions(;
    num_grid_search=250,
    profile_resolution=80,
    min_profile_steps=15,
    num_fine_t=361,
    prediction_interval_resolution=100
);
opt_fnc = (u₀_val, include_perimeter, product, also_compute_separate=false) -> ProfileAnalysisLikelihoodOptions(;
    include_area=true,
    include_perimeter=include_perimeter,
    include_radius=false,
    product=product,
    fix_D=false,
    fix_λ=false,
    fix_u₀=true,
    fixed_D=nothing,
    fixed_λ=nothing,
    fixed_u₀=u₀_val,
    combine=true,
    post=nothing,
    regular_grid_product=regular_grid_product,
    regular_grid_separate=regular_grid_separate,
    also_compute_separate=also_compute_separate
)
wave_opt_fnc = (u₀_val, include_perimeter, product, also_compute_separate=false) -> ProfileAnalysisLikelihoodOptions(;
    include_area=true,
    include_perimeter=include_perimeter,
    include_radius=false,
    product=product,
    fix_D=false,
    fix_λ=false,
    fix_u₀=true,
    fixed_D=nothing,
    fixed_λ=nothing,
    fixed_u₀=u₀_val,
    combine=true,
    post=nothing,
    regular_grid_product=wave_regular_grid_product,
    regular_grid_separate=wave_regular_grid_separate,
    also_compute_separate=also_compute_separate
)

################################################################################################
##
## FULL LIKELIHOOD ANALYSIS
##
################################################################################################

############################################
## Square full likelihood analysis: Area and perimeter
############################################
square_perimeter_area_full_likelihood_results = Vector{Any}(undef, length(u0))
for (idx, u₀_val) in pairs(u0)
    @show idx, u₀_val
    likelihood_options = opt_fnc(u₀_val, true, true, idx == 1) # product argument doesn't matter in this case
    square_perimeter_area_full_likelihood_results[idx] = full_likelihood_analysis(pore_data, save_options, likelihood_options, optimisation_options)
    jldsave(joinpath(jld2_path_full_likelihood, "square_perimeter_area_full_likelihood_results_$idx.jld2"); square_perimeter_area_full_likelihood_results=square_perimeter_area_full_likelihood_results[idx])
end

############################################
## Square full likelihood analysis: Area only
############################################
square_area_full_likelihood_results = Vector{Any}(undef, length(u0))
for (idx, u₀_val) in pairs(u0)
    @show idx, u₀_val
    likelihood_options = opt_fnc(u₀_val, false, true) # product argument doesn't matter in this case
    square_area_full_likelihood_results[idx] = full_likelihood_analysis(pore_data, save_options, likelihood_options, optimisation_options)
    jldsave(joinpath(jld2_path_full_likelihood, "square_area_full_likelihood_results_$idx.jld2"); square_area_full_likelihood_results=square_area_full_likelihood_results[idx])
end

############################################
## Wave full likelihood analysis: Area and perimeter
############################################
wave_perimeter_area_full_likelihood_results = Vector{Any}(undef, length(u0))
for (idx, u₀_val) in pairs(u0)
    @show idx, u₀_val
    likelihood_options = wave_opt_fnc(u₀_val, true, true) # product argument doesn't matter in this case
    wave_perimeter_area_full_likelihood_results[idx] = wave_full_likelihood_analysis(pore_data, save_options, likelihood_options, optimisation_options)
    jldsave(joinpath(jld2_path_full_likelihood, "wave_perimeter_area_full_likelihood_results_$idx.jld2"); wave_perimeter_area_full_likelihood_results=wave_perimeter_area_full_likelihood_results[idx])
end

############################################
## Wave full likelihood analysis: Area only
############################################
wave_area_full_likelihood_results = Vector{Any}(undef, length(u0))
for (idx, u₀_val) in pairs(u0)
    @show idx, u₀_val
    likelihood_options = wave_opt_fnc(u₀_val, false, true) # product argument doesn't matter in this case
    wave_area_full_likelihood_results[idx] = wave_full_likelihood_analysis(pore_data, save_options, likelihood_options, optimisation_options)
    jldsave(joinpath(jld2_path_full_likelihood, "wave_area_full_likelihood_results_$idx.jld2"); wave_area_full_likelihood_results=wave_area_full_likelihood_results[idx])
end

################################################################################################
##
## RUN THE PROFILE ANALYSIS
##
################################################################################################

############################################
## Square profile analysis: Area and perimeter
#############################################
square_perimeter_area_profile_results = Vector{Any}(undef, length(u0))
for (idx, u₀_val) in pairs(u0)
    @show idx, u₀_val
    likelihood_options = opt_fnc(u₀_val, true, true)
    square_perimeter_area_profile_results[idx] = (clear_results ∘ profile_analysis)(pore_data, save_options, likelihood_options, optimisation_options)
end
jldsave(joinpath(jld2_path_profiles, "square_perimeter_area_profile_results.jld2"); square_perimeter_area_results=square_perimeter_area_profile_results)

#############################################
## Square profile analysis: Area only 
#############################################
square_area_profile_results = Vector{Any}(undef, length(u0))
for (idx, u₀_val) in pairs(u0)
    @show idx, u₀_val
    likelihood_options = opt_fnc(u₀_val, false, true)
    square_area_profile_results[idx] = (clear_results ∘ profile_analysis)(pore_data, save_options, likelihood_options, optimisation_options)
end
jldsave(joinpath(jld2_path_profiles, "square_area_results.jld2"); square_area_results=square_area_profile_results)

#############################################
## Wave profile analysis: Area and perimeter
#############################################
wave_perimeter_area_profile_results = Vector{Any}(undef, length(u0))
for (idx, u₀_val) in pairs(u0)
    @show idx, u₀_val
    likelihood_options = opt_fnc(u₀_val, true, true)
    wave_perimeter_area_profile_results[idx] = (clear_results ∘ wave_profile_analysis)(pore_data, save_options, likelihood_options, optimisation_options, geo=:wave)
end
jldsave(joinpath(jld2_path_profiles, "wave_area_perimeter_results.jld2"); wave_area_perimeter_results=wave_perimeter_area_profile_results)

#############################################
## Wave profile analysis: Area only
#############################################
wave_area_profile_results = Vector{Any}(undef, length(u0))
for (idx, u₀_val) in pairs(u0)
    @show idx, u₀_val
    likelihood_options = opt_fnc(u₀_val, idx, false, true)
    wave_area_profile_results[idx] = (clear_results ∘ wave_profile_analysis)(pore_data, save_options, likelihood_options, optimisation_options, geo=:wave)
end
jldsave(joinpath(jld2_path_profiles, "wave_area_results.jld2"); wave_area_results=wave_area_profile_results)

#=
#############################################
## Get results on the square with D and λ separate
#############################################
likelihood_options = opt_fnc(0.2, true, false)
square_area_perimeter_decoupled_profile_results = profile_analysis(pore_data, save_options, likelihood_options, optimisation_options)
cleared_square_area_perimeter_decoupled_profile_results = clear_results(square_area_perimeter_decoupled_profile_results)
jldsave(joinpath(jld2_path_profiles, "square_area_perimeter_decoupled_profile_results.jld2"); cleared_square_area_perimeter_decoupled_profile_results)
=#