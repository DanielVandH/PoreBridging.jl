# Just to make things easier, if I wanted performance I would of course add some types here

"""
    struct ProfileAnalysisLikelihoodOptions

Defines option for the profile likelihood analysis.

# Fields
- `include_area::Bool`: Whether to include the pore area in the likelihood calculation.
- `include_perimeter::Bool`: Whether to include the pore perimeter in the likelihood calculation.
- `include_radius::Bool`: Whether to include the pore radius in the likelihood calculation.
- `product::Bool`: Whether to consider `Dλ` or `(D, λ)`.
- `fix_D::Bool`: Whether to fix `D` in the likelihood calculation.
- `fix_λ::Bool`: Whether to fix `λ` in the likelihood calculation.
- `fix_u₀::Bool`: Whether to fix `u₀` in the likelihood calculation.
- `fixed_D::Union{Nothing, Float64}`: The value of `D` to fix in the likelihood calculation, if `fix_D == true`. Otherwise, `nothing`.
- `fixed_λ::Union{Nothing, Float64}`: The value of `λ` to fix in the likelihood calculation, if `fix_λ == true`. Otherwise, `nothing`.
- `fixed_u₀::Union{Nothing, Float64}`: The value of `u₀` to fix in the likelihood calculation, if `fix_u₀ == true`. Otherwise, `nothing`.
- `combine::Bool`: Whether to combine the experimental data at each time for each summary statistic for computing the noise variances, or to compute an estimate for each time.
- `post::String`: The text to add at the end of the output files. Not used anymore.
- `regular_grid_product`: The regular grid for the `Dλ` parameter space.
- `regular_grid_separate`: The regular grid for the `(D, λ)` parameter space.
"""
Base.@kwdef struct ProfileAnalysisLikelihoodOptions
    include_area
    include_perimeter
    include_radius
    product
    fix_D
    fix_λ
    fix_u₀
    fixed_D
    fixed_λ
    fixed_u₀
    combine
    post
    regular_grid_product
    regular_grid_separate
    also_compute_separate
end

"""
    struct ProfileAnalysisSaveOptions

Defines options for saving the results of the profile likelihood analysis.

# Fields
- `save_path::String`: The path to save the results to.
- `initial_time::Real`: The initial time of the simulation.
- `final_time::Real`: The final time of the simulation.
- `saveat_times::Vector{<:Real}`: The times to save the numerical solutions at for the objective function.
- `snapshot_save_times::Vector{<:Real}`: The times to save the snapshots at.
"""
Base.@kwdef struct ProfileAnalysisSaveOptions
    save_path
    initial_time
    final_time
    saveat_times
    snapshot_save_times
end

"""
    struct PoreDataOptions 

Defines options for the pore data.

# Fields
- `data`: The data set.
- `mesh`: The square mesh.
- `pore_area`: The square pore area.
- `pore_perimeter`: The square pore perimeter.
- `pore_smallest_rad`: The smallest distance from the `centroid` to the square pore boundary.
- `centroid`: The centroid of the square pore.
- `wave_mesh`: The mesh for the wave pore.
- `wave_pore_area`: The wave pore area.
- `wave_pore_perimeter`: The wave pore perimeter.
- `wave_pore_smallest_rad`: The smallest distance from the `wave_centroid` to the wave pore boundary.
- `wave_centroid`: The centroid of the wave pore.
"""
Base.@kwdef struct PoreDataOptions
    data
    mesh
    pore_area
    pore_perimeter
    pore_smallest_rad
    centroid
    wave_mesh
    wave_pore_area
    wave_pore_perimeter
    wave_pore_smallest_rad
    wave_centroid
end

"""
    struct ProfileAnalysisResults

Defines the results of the profile likelihood analysis when considering the 
square geometry.

# Fields
- `Clikprob`: The likelihood problem.
- `Cgs_ir`: The grid search results. (The `ir` means irregular grid.)
- `Cloglik_vals`: The log-likelihood values from the grid search.
- `Crefined_likprob`: The refined likelihood problem, using the results from the grid search. 
- `Cmle_sol`: The maximum likelihood solution.
- `Cprof`: The profile likelihood solution.
- `Cfig_profiles`: The figure of the profiles.
- `Cfig_snapshots`: The figure of the snapshots.
- `C_D`: The `D` values for the profile likelihood analysis used for the snapshots.
- `C_λ`: The `λ` values for the profile likelihood analysis used for the snapshots.
- `C_u₀`: The `u₀` values for the profile likelihood analysis used for the snapshots.
- `Cindividual_intervals`: The individual intervals for each parameter for the prediction intervals.
- `Cunion_intervals`: The unions of the individual intervals for each parameter for the prediction intervals.
- `Cq_vals`: The `q` values for the prediction intervals, meaning the prediction function values.
- `Cparam_ranges`: The parameter ranges for the prediction intervals.
- `Cq_mle`: The value of the prediction function at the MLE.
- `Csaveat_times_fine`: The times to save the numerical solutions at for the objective function when considering the prediction intervals. 
- `Cfig_predictions`: The figure of the prediction intervals.
- `Cindiv_idx`: The indices of the individual intervals for each parameter for the prediction intervals.
- `Cfig_closing_times`: The figure of the closing times.
- `Cfig_wave_snapshots`: The figure of the wave snapshots.
- `Cwave_individual_intervals`: The individual intervals for each parameter for the wave prediction intervals.
- `Cwave_union_intervals`: The unions of the individual intervals for each parameter for the wave prediction intervals.
- `Cwave_q_vals`: The `q` values for the wave prediction intervals, meaning the prediction function values.
- `Cwave_param_ranges`: The parameter ranges for the wave prediction intervals.
- `Cwave_q_mle`: The value of the prediction function at the MLE for the wave prediction intervals.
- `Cwave_saveat_times_fine`: The times to save the numerical solutions at for the objective function when considering the wave prediction intervals.
- `Cfig_wave_predictions`: The figure of the wave prediction intervals.
- `Cwave_indiv_idx`: The indices of the individual intervals for each parameter for the wave prediction intervals.
- `Cfig_wave_closing_times`: The figure of the wave closing times.
- `pore_data`: The pore data; see [`PoreDataOptions`](@ref).
- `save_options`: The save options; see [`ProfileAnalysisSaveOptions`](@ref).
- `likelihood_options`: The likelihood options; see [`ProfileAnalysisLikelihoodOptions`](@ref).
- `optimisation_options`: The optimisation options; see [`OptimisationOptions`](@ref).
"""
Base.@kwdef struct ProfileAnalysisResults
    Clikprob
    Cgs_ir
    Cloglik_vals
    Crefined_likprob
    Cmle_sol
    Cprof
    Cfig_profiles
    Cfig_snapshots
    C_D
    C_λ
    C_u₀
    Cindividual_intervals
    Cunion_intervals
    Cq_vals
    Cparam_ranges
    Cq_mle
    Csaveat_times_fine
    Cfig_predictions
    Cindiv_idx
    Cfig_closing_times
    Cfig_wave_snapshots
    Cwave_individual_intervals
    Cwave_union_intervals
    Cwave_q_vals
    Cwave_param_ranges
    Cwave_q_mle
    Cwave_saveat_times_fine
    Cfig_wave_predictions
    Cwave_indiv_idx
    Cfig_wave_closing_times
    pore_data
    save_options
    likelihood_options
    optimisation_options
end

"""
    struct FullLikelihoodResults

Defines the results of the profile likelihood analysis when considering the 
square geometry.

# Fields
- `likprob_product`: The likelihood problem for the `Dλ` parameter space.
- `likprob_separate`: The likelihood problem for the `(D, λ)` parameter space.
- `refined_likprob_product`: The refined likelihood problem for the `Dλ` parameter space.
- `refined_likprob_separate`: The refined likelihood problem for the `(D, λ)` parameter space.
- `loglik_vals_product`: The normalised log-likelihood values for the `Dλ` parameter space.
- `loglik_vals_separate`: The normalised log-likelihood values for the `(D, λ)` parameter space.
- `mle_sol_product`: The MLE for the `Dλ` parameter space.
- `mle_sol_separate`: The MLE for the `(D, λ)` parameter space.
- `intervals`: The prediction intervals.
- `q_vals`: The `q` values for the prediction intervals, meaning the prediction function values.
- `param_ranges`: The parameter ranges for the prediction intervals.
- `q_mle`: The value of the prediction function at the MLE.
- `saveat_times_fine`: The times to save the numerical solutions at for the objective function when considering the prediction intervals.
- `wave_intervals`: The prediction intervals for the wave geometry.
- `wave_q_vals`: The `q` values for the wave prediction intervals, meaning the prediction function values.
- `wave_param_ranges`: The parameter ranges for the wave prediction intervals.
- `wave_q_mle`: The value of the prediction function at the MLE for the wave prediction intervals.
- `wave_saveat_times_fine`: The times to save the numerical solutions at for the objective function when considering the wave prediction intervals.
- `pore_data`: The pore data; see [`PoreDataOptions`](@ref).
- `save_options`: The save options; see [`ProfileAnalysisSaveOptions`](@ref).
- `likelihood_options`: The likelihood options; see [`ProfileAnalysisLikelihoodOptions`](@ref).
- `optimisation_options`: The optimisation options; see [`OptimisationOptions`](@ref).
"""
Base.@kwdef struct FullLikelihoodResults
    likprob_product
    likprob_separate
    refined_likprob_product
    refined_likprob_separate
    loglik_vals_product
    loglik_vals_separate
    mle_sol_product
    mle_sol_separate
    intervals
    q_vals
    param_ranges
    q_mle
    saveat_times_fine
    wave_intervals
    wave_q_vals
    wave_param_ranges
    wave_q_mle
    wave_saveat_times_fine
    pore_data
    save_options
    likelihood_options
    optimisation_options
end

"""
    struct ReverseProfileAnalysisResults

Defines the results of the profile likelihood analysis when considering the 
wave geometry. See [`ProfileAnalysisResults`](@ref), but just interchange `square` 
and `wave`.
"""
Base.@kwdef struct ReverseProfileAnalysisResults
    Clikprob
    Cgs_ir
    Cloglik_vals
    Crefined_likprob
    Cmle_sol
    Cprof
    Cfig_profiles
    Cfig_snapshots
    C_D
    C_λ
    C_u₀
    Cindividual_intervals
    Cunion_intervals
    Cq_vals
    Cparam_ranges
    Cq_mle
    Csaveat_times_fine
    Cfig_predictions
    Cindiv_idx
    Cfig_closing_times
    Cfig_square_snapshots
    Csquare_individual_intervals
    Csquare_union_intervals
    Csquare_q_vals
    Csquare_param_ranges
    Csquare_q_mle
    Csquare_saveat_times_fine
    Cfig_square_predictions
    Csquare_indiv_idx
    Cfig_square_closing_times
    pore_data
    save_options
    likelihood_options
    optimisation_options
end

"""
    struct ReverseFullLikelihoodResults

Defines the results of the profile likelihood analysis when considering the
wave geometry. See [`FullLikelihoodResults`](@ref), but just interchange `square`
and `wave`.
"""
Base.@kwdef struct ReverseFullLikelihoodResults
    likprob_product
    likprob_separate
    refined_likprob_product
    refined_likprob_separate
    loglik_vals_product
    loglik_vals_separate
    mle_sol_product
    mle_sol_separate
    intervals
    q_vals
    param_ranges
    q_mle
    saveat_times_fine
    square_intervals
    square_q_vals
    square_param_ranges
    square_q_mle
    square_saveat_times_fine
    pore_data
    save_options
    likelihood_options
    optimisation_options
end

"""
    struct OptimisationOptions

Defines options for the optimisation.

# Fields
- `num_grid_search::Int`: The number of grid search points to use.
- `profile_resolution::Int`: The number of points to use in each direction when computing the profile likelihoods. 
- `min_profile_steps::Int`: The minimum number of steps to use when computing the profile likelihoods.
- `num_fine_t::Int`: The number of points to use when computing the prediction intervals.
- `prediction_interval_resolution::Int`: The number of samples to compute for each parameter when computing the prediction intervals.
"""
Base.@kwdef struct OptimisationOptions
    num_grid_search
    profile_resolution
    min_profile_steps
    num_fine_t
    prediction_interval_resolution
end