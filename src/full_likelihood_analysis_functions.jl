"""
    full_likelihood_analysis(pore_data, save_options, likelihood_options, optimisation_options)

Performs the full likelihood analysis on the square geometry.

# Arguments
- `pore_data`: The [`PoreDataOptions`](@ref).
- `save_options`: The [`ProfileAnalysisSaveOptions](@ref)`.
- `likelihood_options`: The [`ProfileLikelihoodOptions](@ref)`.
- `optimisation_options`: The [`ProfileOptimisationOptions](@ref)`.

# Outputs
- `results`: The [`FullLikelihoodResults`](@ref).
"""
function full_likelihood_analysis(pore_data, save_options, likelihood_options, optimisation_options)
    (;
        data,
        mesh,
        pore_area,
        pore_perimeter,
        pore_smallest_rad,
        centroid,
        wave_mesh,
        wave_pore_area,
        wave_pore_perimeter,
        wave_pore_smallest_rad,
        wave_centroid) = pore_data
    (;
        save_path,
        initial_time,
        final_time,
        saveat_times,
        snapshot_save_times) = save_options
    (;
        include_area,
        include_perimeter,
        include_radius,
        product,
        fix_D,
        fix_λ,
        fix_u₀,
        fixed_u₀,
        fixed_λ,
        fixed_D,
        combine,
        post,
        regular_grid_product,
        regular_grid_separate,
        also_compute_separate) = likelihood_options
    (;
        num_grid_search,
        profile_resolution,
        min_profile_steps,
        num_fine_t,
        prediction_interval_resolution) = optimisation_options

    areas, perimeters, radii, σ = get_pore_data(data, saveat_times; combine=combine)
    wave_areas, wave_perimeters, wave_radii, wave_σ = get_pore_data(data, saveat_times, "SWV500"; combine=true)

    likprob_product = construct_likelihood_problem(
        mesh,
        initial_time,
        final_time,
        saveat_times,
        areas,
        perimeters,
        radii,
        σ,
        pore_area,
        pore_perimeter,
        pore_smallest_rad,
        centroid;
        include_area,
        include_perimeter,
        include_radius,
        product=true,
        fix_D,
        fix_λ,
        fix_u₀,
        fixed_D,
        fixed_u₀,
        fixed_λ,
        geo=:square
    )
    likprob_separate = construct_likelihood_problem(
        mesh,
        initial_time,
        final_time,
        saveat_times,
        areas,
        perimeters,
        radii,
        σ,
        pore_area,
        pore_perimeter,
        pore_smallest_rad,
        centroid;
        include_area,
        include_perimeter,
        include_radius,
        product=false,
        fix_D,
        fix_λ,
        fix_u₀,
        fixed_D,
        fixed_u₀,
        fixed_λ,
        geo=:square
    )

    @time gs_product, loglik_vals_product = perform_grid_search(likprob_product, regular_grid_product)
    if also_compute_separate
        @time gs_separate, loglik_vals_separate = perform_grid_search(likprob_separate, regular_grid_separate)
    else
        gs_separate = nothing
        loglik_vals_separate = nothing
    end

    @time refined_likprob_product, mle_product = find_mle(likprob_product, gs_product)
    if also_compute_separate
        @time refined_likprob_separate, mle_separate = find_mle(likprob_separate, gs_separate)
    else
        refined_likprob_separate = nothing
        mle_separate = nothing
    end

    normalised_loglik_vals_product = loglik_vals_product .- get_maximum(mle_product)
    if also_compute_separate
        normalised_loglik_vals_separate = loglik_vals_separate .- get_maximum(mle_separate)
    else
        normalised_loglik_vals_separate = nothing
    end
    feasible_idx = findall(>(ProfileLikelihood.get_chisq_threshold(0.95, 2)), normalised_loglik_vals_product)
    grid_1 = get_range(regular_grid_product, 1)
    grid_2 = get_range(regular_grid_product, 2)
    feasible_pts = [(grid_1[I[1]], grid_2[I[2]]) for I in feasible_idx]

    nt = num_fine_t
    full_q_vals = zeros(3nt + 1, length(feasible_pts))
    parallel = true
    ext_initial_time = initial_time
    ext_final_time = 5final_time
    saveat_times_fine = LinRange(ext_initial_time, ext_final_time, nt)
    p = (mesh=mesh,
        ext_initial_time=ext_initial_time,
        ext_final_time=ext_final_time,
        saveat_times_fine=saveat_times_fine,
        pore_area=pore_area,
        pore_perimeter=pore_perimeter,
        pore_smallest_rad=pore_smallest_rad,
        centroid=centroid,
        use_concave_hull=false)
    fnc = construct_prediction_function(refined_likprob_product)
    Base.Threads.@threads for (i, θ) in collect(enumerate(feasible_pts))
        @show i, θ
        full_q_vals[:, i] .= fnc(θ, p)
    end
    q_mle = fnc(get_mle(mle_product), p)
    q_lwr = minimum(full_q_vals; dims=2) |> vec
    q_upr = maximum(full_q_vals; dims=2) |> vec
    intervals = [(l, u) for (l, u) in zip(q_lwr, q_upr)]

    wave_full_q_vals = zeros(3nt + 1, length(feasible_pts))
    wave_saveat_times_fine = LinRange(ext_initial_time, ext_final_time, nt)
    p = (mesh=wave_mesh,
        ext_initial_time=ext_initial_time,
        ext_final_time=ext_final_time,
        saveat_times_fine=wave_saveat_times_fine,
        pore_area=wave_pore_area,
        pore_perimeter=wave_pore_perimeter,
        pore_smallest_rad=wave_pore_smallest_rad,
        centroid=wave_centroid,
        use_concave_hull=false)
    fnc = construct_prediction_function(refined_likprob_product)
    Base.Threads.@threads for (i, θ) in collect(enumerate(feasible_pts))
        @show i, θ
        wave_full_q_vals[:, i] .= fnc(θ, p)
    end
    wave_q_mle = fnc(get_mle(mle_product), p)
    wave_q_lwr = minimum(wave_full_q_vals; dims=2) |> vec
    wave_q_upr = maximum(wave_full_q_vals; dims=2) |> vec
    wave_intervals = [(l, u) for (l, u) in zip(wave_q_lwr, wave_q_upr)]

    return FullLikelihoodResults(;
        likprob_product=likprob_product,
        likprob_separate=likprob_separate,
        refined_likprob_product=refined_likprob_product,
        refined_likprob_separate=refined_likprob_separate,
        loglik_vals_product=normalised_loglik_vals_product,
        loglik_vals_separate=normalised_loglik_vals_separate,
        mle_sol_product=mle_product,
        mle_sol_separate=mle_separate,
        intervals=intervals,
        q_vals=full_q_vals,
        param_ranges=(grid_1, grid_2),
        q_mle=q_mle,
        saveat_times_fine=saveat_times_fine,
        wave_intervals=wave_intervals,
        wave_q_vals=wave_full_q_vals,
        wave_param_ranges=(grid_1, grid_2),
        wave_q_mle=wave_q_mle,
        wave_saveat_times_fine=wave_saveat_times_fine,
        pore_data,
        save_options,
        likelihood_options,
        optimisation_options
    )
end

"""
    wave_full_likelihood_analysis(pore_data, save_options, likelihood_options, optimisation_options)

Performs the full likelihood analysis on the wave geometry.

# Arguments
- `pore_data`: The [`PoreDataOptions`](@ref).
- `save_options`: The [`ProfileAnalysisSaveOptions](@ref)`.
- `likelihood_options`: The [`ProfileLikelihoodOptions](@ref)`.
- `optimisation_options`: The [`ProfileOptimisationOptions](@ref)`.

# Outputs
- `results`: The [`FullLikelihoodResults`](@ref).
"""
function wave_full_likelihood_analysis(pore_data, save_options, likelihood_options, optimisation_options)
    (;
        data,
        wave_mesh,
        wave_pore_area,
        wave_pore_perimeter,
        wave_pore_smallest_rad,
        wave_centroid) = pore_data
    (;
        save_path,
        initial_time,
        final_time,
        saveat_times,
        snapshot_save_times) = save_options
    (;
        include_area,
        include_perimeter,
        include_radius,
        product,
        fix_D,
        fix_λ,
        fix_u₀,
        fixed_u₀,
        fixed_λ,
        fixed_D,
        combine,
        post,
        regular_grid_product,
        regular_grid_separate,
        also_compute_separate) = likelihood_options
    (;
        num_grid_search,
        profile_resolution,
        min_profile_steps,
        num_fine_t,
        prediction_interval_resolution) = optimisation_options

    areas, perimeters, radii, σ = get_pore_data(data, saveat_times; combine=combine)
    wave_areas, wave_perimeters, wave_radii, wave_σ = get_pore_data(data, saveat_times, "SWV500"; combine=true)

    likprob_product = construct_likelihood_problem(
        wave_mesh,
        initial_time,
        final_time,
        saveat_times,
        wave_areas,
        wave_perimeters,
        wave_radii,
        wave_σ,
        wave_pore_area,
        wave_pore_perimeter,
        wave_pore_smallest_rad,
        wave_centroid;
        include_area,
        include_perimeter,
        include_radius,
        product=true,
        fix_D,
        fix_λ,
        fix_u₀,
        fixed_D,
        fixed_u₀,
        fixed_λ,
        geo=:wave
    )
    likprob_separate = construct_likelihood_problem(
        wave_mesh,
        initial_time,
        final_time,
        saveat_times,
        wave_areas,
        wave_perimeters,
        wave_radii,
        wave_σ,
        wave_pore_area,
        wave_pore_perimeter,
        wave_pore_smallest_rad,
        wave_centroid;
        include_area,
        include_perimeter,
        include_radius,
        product=false,
        fix_D,
        fix_λ,
        fix_u₀,
        fixed_D,
        fixed_u₀,
        fixed_λ,
        geo=:wave
    )

    @time gs_product, loglik_vals_product = perform_grid_search(likprob_product, regular_grid_product)
    if also_compute_separate
        @time gs_separate, loglik_vals_separate = perform_grid_search(likprob_separate, regular_grid_separate)
    else
        gs_separate = nothing
        loglik_vals_separate = nothing
    end

    @time refined_likprob_product, mle_product = find_mle(likprob_product, gs_product)
    if also_compute_separate
        @time refined_likprob_separate, mle_separate = find_mle(likprob_separate, gs_separate)
    else
        refined_likprob_separate = nothing
        mle_separate = nothing
    end

    normalised_loglik_vals_product = loglik_vals_product .- get_maximum(mle_product)
    if also_compute_separate
        normalised_loglik_vals_separate = loglik_vals_separate .- get_maximum(mle_separate)
    else
        normalised_loglik_vals_separate = nothing
    end
    feasible_idx = findall(>(ProfileLikelihood.get_chisq_threshold(0.95, 2)), normalised_loglik_vals_product)
    grid_1 = get_range(regular_grid_product, 1)
    grid_2 = get_range(regular_grid_product, 2)
    feasible_pts = [(grid_1[I[1]], grid_2[I[2]]) for I in feasible_idx]
    @show length(feasible_pts)

    nt = num_fine_t
    full_q_vals = zeros(3nt + 1, length(feasible_pts))
    parallel = true
    ext_initial_time = initial_time
    ext_final_time = 5final_time
    saveat_times_fine = LinRange(ext_initial_time, ext_final_time, nt)
    p = (mesh=wave_mesh,
        ext_initial_time=ext_initial_time,
        ext_final_time=ext_final_time,
        saveat_times_fine=saveat_times_fine,
        pore_area=wave_pore_area,
        pore_perimeter=wave_pore_perimeter,
        pore_smallest_rad=wave_pore_smallest_rad,
        centroid=wave_centroid,
        use_concave_hull=false)
    fnc = construct_prediction_function(refined_likprob_product)
    Base.Threads.@threads for (i, θ) in collect(enumerate(feasible_pts))
        @show i, θ
        full_q_vals[:, i] .= fnc(θ, p)
    end
    q_mle = fnc(get_mle(mle_product), p)
    q_lwr = minimum(full_q_vals; dims=2) |> vec
    q_upr = maximum(full_q_vals; dims=2) |> vec
    intervals = [(l, u) for (l, u) in zip(q_lwr, q_upr)]

    (;
        mesh,
        pore_area,
        pore_perimeter,
        pore_smallest_rad,
        centroid) = pore_data
    areas, perimeters, radii, σ = get_pore_data(data, saveat_times; combine=combine)

    square_full_q_vals = zeros(3nt + 1, length(feasible_pts))
    square_saveat_times_fine = LinRange(ext_initial_time, ext_final_time, nt)
    p = (mesh=mesh,
        ext_initial_time=ext_initial_time,
        ext_final_time=ext_final_time,
        saveat_times_fine=square_saveat_times_fine,
        pore_area=pore_area,
        pore_perimeter=pore_perimeter,
        pore_smallest_rad=pore_smallest_rad,
        centroid=centroid,
        use_concave_hull=false)
    fnc = construct_prediction_function(refined_likprob_product)
    Base.Threads.@threads for (i, θ) in collect(enumerate(feasible_pts))
        @show i, θ
        square_full_q_vals[:, i] .= fnc(θ, p)
    end
    square_q_mle = fnc(get_mle(mle_product), p)
    square_q_lwr = minimum(square_full_q_vals; dims=2) |> vec
    square_q_upr = maximum(square_full_q_vals; dims=2) |> vec
    square_intervals = [(l, u) for (l, u) in zip(square_q_lwr, square_q_upr)]

    return ReverseFullLikelihoodResults(;
        likprob_product=likprob_product,
        likprob_separate=likprob_separate,
        refined_likprob_product=refined_likprob_product,
        refined_likprob_separate=refined_likprob_separate,
        loglik_vals_product=normalised_loglik_vals_product,
        loglik_vals_separate=normalised_loglik_vals_separate,
        mle_sol_product=mle_product,
        mle_sol_separate=mle_separate,
        intervals=intervals,
        q_vals=full_q_vals,
        param_ranges=(grid_1, grid_2),
        q_mle=q_mle,
        saveat_times_fine=saveat_times_fine,
        square_intervals=square_intervals,
        square_q_vals=square_full_q_vals,
        square_param_ranges=(grid_1, grid_2),
        square_q_mle=square_q_mle,
        square_saveat_times_fine=square_saveat_times_fine,
        pore_data,
        save_options,
        likelihood_options,
        optimisation_options
    )
end