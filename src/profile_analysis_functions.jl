"""
    profile_analysis(pore_data, save_options, likelihood_options, optimisation_options)

Performs the profile likelihood analysis on the square geometry.

# Arguments
- `pore_data`: The [`PoreDataOptions`](@ref).
- `save_options`: The [`ProfileAnalysisSaveOptions](@ref)`.
- `likelihood_options`: The [`ProfileLikelihoodOptions](@ref)`.
- `optimisation_options`: The [`ProfileOptimisationOptions](@ref)`.

# Outputs
- `results`: The [`ProfileAnalysisResults`](@ref).
"""
function profile_analysis(pore_data, save_options, likelihood_options, optimisation_options)
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
        post) = likelihood_options
    (;
        num_grid_search,
        profile_resolution,
        min_profile_steps,
        num_fine_t,
        prediction_interval_resolution) = optimisation_options

    areas, perimeters, radii, σ = get_pore_data(data, saveat_times; combine=combine)
    wave_areas, wave_perimeters, wave_radii, wave_σ = get_pore_data(data, saveat_times, "SWV500"; combine=true)

    likprob = construct_likelihood_problem(
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
        product,
        fix_D,
        fix_λ,
        fix_u₀,
        fixed_D,
        fixed_u₀,
        fixed_λ,
        geo=:square
    )

    gs_ir, loglik_vals = perform_grid_search(likprob, num_grid_search)
    refined_likprob, mle_sol = find_mle(likprob, gs_ir)

    prof = find_profiles(refined_likprob, mle_sol, profile_resolution, min_profile_steps)
    fig_profiles = plot_profiles(prof; nrow=1, ncol=ProfileLikelihood.number_of_parameters(likprob),
        latex_names=[L"%$s" for s in ProfileLikelihood.get_syms(likprob)],
        fig_kwargs=(fontsize=38, resolution=(2158.8823f0, 470.17322f0)),
        axis_kwargs=(width=600, height=300),
        show_points=true,
        spline=true)
    resize_to_layout!(fig_profiles)

    fig_snapshots, _D, _λ, _u₀ = plot_snapshots(prof, snapshot_save_times;
        product, fix_D, fix_λ, fix_u₀, fixed_D, fixed_λ, fixed_u₀)
    resize_to_layout!(fig_snapshots)

    individual_intervals, union_intervals, q_vals, param_ranges, q_mle, saveat_times_fine =
        compute_prediction_intervals(prof, num_fine_t; resolution=prediction_interval_resolution)
    fig_predictions, indiv_idx = plot_prediction_functions(saveat_times_fine,
        individual_intervals,
        union_intervals,
        q_mle,
        areas,
        perimeters,
        radii,
        saveat_times,
        product,
        fix_D,
        fix_λ,
        fix_u₀)
    resize_to_layout!(fig_predictions)

    fig_closing_times = plot_closing_time_densities(q_vals, q_mle, indiv_idx)
    resize_to_layout!(fig_closing_times)

    fig_wave_snapshots = plot_snapshots(_D, _λ, _u₀, wave_mesh, snapshot_save_times, wave_centroid;
        use_concave_hull=true,
        rotate=true,
        shift_vals=true)
    resize_to_layout!(fig_wave_snapshots)

    wave_individual_intervals, wave_union_intervals, wave_q_vals, wave_param_ranges, wave_q_mle, wave_saveat_times_fine =
        compute_prediction_intervals(prof, num_fine_t;
            mesh=wave_mesh,
            pore_area=wave_pore_area,
            pore_perimeter=wave_pore_perimeter,
            pore_smallest_rad=wave_pore_smallest_rad,
            centroid=wave_centroid,
            use_concave_hull=false,
            resolution=prediction_interval_resolution)
    fig_wave_predictions, wave_indiv_idx = plot_prediction_functions(wave_saveat_times_fine,
        wave_individual_intervals,
        wave_union_intervals,
        wave_q_mle,
        wave_areas,
        wave_perimeters,
        wave_radii,
        saveat_times,
        product,
        fix_D,
        fix_λ,
        fix_u₀)
    resize_to_layout!(fig_wave_predictions)

    fig_wave_closing_times = plot_closing_time_densities(wave_q_vals, wave_q_mle, wave_indiv_idx)
    resize_to_layout!(fig_wave_closing_times)

    results = ProfileAnalysisResults(;
        Clikprob=likprob,
        Cgs_ir=gs_ir,
        Cloglik_vals=loglik_vals,
        Crefined_likprob=refined_likprob,
        Cmle_sol=mle_sol,
        Cprof=prof,
        Cfig_profiles=fig_profiles,
        Cfig_snapshots=fig_snapshots,
        C_D=_D,
        C_λ=_λ,
        C_u₀=_u₀,
        Cindividual_intervals=individual_intervals,
        Cunion_intervals=union_intervals,
        Cq_vals=q_vals,
        Cparam_ranges=param_ranges,
        Cq_mle=q_mle,
        Csaveat_times_fine=saveat_times_fine,
        Cfig_predictions=fig_predictions,
        Cindiv_idx=indiv_idx,
        Cfig_closing_times=fig_closing_times,
        Cfig_wave_snapshots=fig_wave_snapshots,
        Cwave_individual_intervals=wave_individual_intervals,
        Cwave_union_intervals=wave_union_intervals,
        Cwave_q_vals=wave_q_vals,
        Cwave_param_ranges=wave_param_ranges,
        Cwave_q_mle=wave_q_mle,
        Cwave_saveat_times_fine=wave_saveat_times_fine,
        Cfig_wave_predictions=fig_wave_predictions,
        Cwave_indiv_idx=wave_indiv_idx,
        Cfig_wave_closing_times=fig_wave_closing_times,
        pore_data=pore_data,
        save_options=save_options,
        likelihood_options=likelihood_options,
        optimisation_options=optimisation_options
    )
    return results
end

"""
    wave_profile_analysis(pore_data, save_options, likelihood_options, optimisation_options)

Performs the profile likelihood analysis on the wave geometry.

# Arguments
- `pore_data`: The [`PoreDataOptions`](@ref).
- `save_options`: The [`ProfileAnalysisSaveOptions](@ref)`.
- `likelihood_options`: The [`ProfileLikelihoodOptions](@ref)`.
- `optimisation_options`: The [`ProfileOptimisationOptions](@ref)`.

# Outputs
- `results`: The [`ReverseProfileAnalysisResults`](@ref).
"""
function wave_profile_analysis(pore_data, save_options, likelihood_options, optimisation_options)
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
        post) = likelihood_options
    (;
        num_grid_search,
        profile_resolution,
        min_profile_steps,
        num_fine_t,
        prediction_interval_resolution) = optimisation_options

    areas, perimeters, radii, σ = get_pore_data(data, saveat_times; combine=combine)
    wave_areas, wave_perimeters, wave_radii, wave_σ = get_pore_data(data, saveat_times, "SWV500"; combine=true)

    likprob = construct_likelihood_problem(
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
        product,
        fix_D,
        fix_λ,
        fix_u₀,
        fixed_D,
        fixed_u₀,
        fixed_λ,
        geo=:wave
    )

    gs_ir, loglik_vals = perform_grid_search(likprob, num_grid_search)
    refined_likprob, mle_sol = find_mle(likprob, gs_ir)

    prof = find_profiles(refined_likprob, mle_sol, profile_resolution, min_profile_steps)
    fig_profiles = plot_profiles(prof; nrow=1, ncol=ProfileLikelihood.number_of_parameters(likprob),
        latex_names=[L"%$s" for s in ProfileLikelihood.get_syms(likprob)],
        fig_kwargs=(fontsize=38, resolution=(2158.8823f0, 470.17322f0)),
        axis_kwargs=(width=600, height=300),
        show_points=true,
        spline=true)
    resize_to_layout!(fig_profiles)

    fig_snapshots, _D, _λ, _u₀ = plot_snapshots(prof, snapshot_save_times;
        product, fix_D, fix_λ, fix_u₀, fixed_D, fixed_λ, fixed_u₀,
        use_concave_hull=true, rotate=true, shift_vals=true)
    resize_to_layout!(fig_snapshots)

    individual_intervals, union_intervals, q_vals, param_ranges, q_mle, saveat_times_fine =
        compute_prediction_intervals(prof, num_fine_t; resolution=prediction_interval_resolution)
    fig_predictions, indiv_idx = plot_prediction_functions(saveat_times_fine,
        individual_intervals,
        union_intervals,
        q_mle,
        wave_areas,
        wave_perimeters,
        wave_radii,
        saveat_times,
        product,
        fix_D,
        fix_λ,
        fix_u₀)
    resize_to_layout!(fig_predictions)

    fig_closing_times = plot_closing_time_densities(q_vals, q_mle, indiv_idx)
    resize_to_layout!(fig_closing_times)

    fig_square_snapshots = plot_snapshots(_D, _λ, _u₀, mesh, snapshot_save_times, centroid;
        use_concave_hull=false,
        rotate=false,
        shift_vals=false)
    resize_to_layout!(fig_square_snapshots)

    square_individual_intervals, square_union_intervals, square_q_vals, square_param_ranges, square_q_mle, square_saveat_times_fine =
        compute_prediction_intervals(prof, num_fine_t;
            mesh=mesh,
            pore_area=pore_area,
            pore_perimeter=pore_perimeter,
            pore_smallest_rad=pore_smallest_rad,
            centroid=centroid,
            use_concave_hull=false,
            resolution=prediction_interval_resolution)
    fig_square_predictions, square_indiv_idx = plot_prediction_functions(square_saveat_times_fine,
        square_individual_intervals,
        square_union_intervals,
        square_q_mle,
        areas,
        perimeters,
        radii,
        saveat_times,
        product,
        fix_D,
        fix_λ,
        fix_u₀)
    resize_to_layout!(fig_square_predictions)

    fig_square_closing_times = plot_closing_time_densities(square_q_vals, square_q_mle, square_indiv_idx)
    resize_to_layout!(fig_square_closing_times)

    results = ReverseProfileAnalysisResults(;
        Clikprob=likprob,
        Cgs_ir=gs_ir,
        Cloglik_vals=loglik_vals,
        Crefined_likprob=refined_likprob,
        Cmle_sol=mle_sol,
        Cprof=prof,
        Cfig_profiles=fig_profiles,
        Cfig_snapshots=fig_snapshots,
        C_D=_D,
        C_λ=_λ,
        C_u₀=_u₀,
        Cindividual_intervals=individual_intervals,
        Cunion_intervals=union_intervals,
        Cq_vals=q_vals,
        Cparam_ranges=param_ranges,
        Cq_mle=q_mle,
        Csaveat_times_fine=saveat_times_fine,
        Cfig_predictions=fig_predictions,
        Cindiv_idx=indiv_idx,
        Cfig_closing_times=fig_closing_times,
        Cfig_square_snapshots=fig_square_snapshots,
        Csquare_individual_intervals=square_individual_intervals,
        Csquare_union_intervals=square_union_intervals,
        Csquare_q_vals=square_q_vals,
        Csquare_param_ranges=square_param_ranges,
        Csquare_q_mle=square_q_mle,
        Csquare_saveat_times_fine=square_saveat_times_fine,
        Cfig_square_predictions=fig_square_predictions,
        Csquare_indiv_idx=square_indiv_idx,
        Cfig_square_closing_times=fig_square_closing_times,
        pore_data=pore_data,
        save_options=save_options,
        likelihood_options=likelihood_options,
        optimisation_options=optimisation_options
    )
    return results
end


function plot_square_likelihood_analysis(all_results, u0, pore_data, save_options, profile_results;
    square_ptb_xticks=(18:4:34, [L"%$s" for s in 18:4:34]),
    wave_ptb_xticks=(11:4:28, [L"%$s" for s in 11:4:28]))
    colors = [:black, :blue, :magenta]
    levels = [ProfileLikelihood.get_chisq_threshold(0.95, 2)]
    threshold = ProfileLikelihood.get_chisq_threshold(0.95, 1)
    linewidth = 4
    markersize = 27
    data_markercolor = :blue
    strokearound = false
    strokewidth = 4
    areas, perimeters, radii, σ = get_pore_data(pore_data.data, save_options.saveat_times; combine=true)
    wave_areas, wave_perimeters, wave_radii, wave_σ = get_pore_data(pore_data.data, save_options.saveat_times[2], "SWV500"; combine=true)
    days = [repeat([save_options.saveat_times[i]], length(areas[i])) for i in eachindex(save_options.saveat_times)]
    wave_days = repeat([save_options.saveat_times[2]], length(wave_areas))

    fig1 = Figure(fontsize=63)
    ax1 = Axis(fig1[1, 1],
        xlabel=L"$D\lambda$ [$\mu$m$^2$/day$^2$]",
        ylabel=L"$\lambda$ [1/day]",
        title=L"(a):$ $ Likelihood CRs",
        titlealign=:left,
        xticks=(0:100:400, [L"%$s" for s in 0:100:400]),
        yticks=(0:2.5:5, [L"%$s" for s in 0:2.5:5.0]),
        width=600,
        height=400)
    for (i, result) in enumerate(all_results)
        grid_1 = result.param_ranges[1]
        grid_2 = result.param_ranges[2]
        normalised_loglik_vals_product = result.loglik_vals_product
        contour!(ax1, grid_1, grid_2, normalised_loglik_vals_product, levels=levels, color=colors[i], linewidth=linewidth)
        D̂ = profile_results[i].mle_val
        vlines!(ax1, [D̂], linestyle=:dash, linewidth=linewidth, color=colors[i])
    end
    xlims!(ax1, 0, 400)

    ax2 = Axis(fig1[2, 1],
        xlabel=L"$D\lambda$ [$\mu$m$^2$/day$^2$]",
        ylabel=L"\ell_p(D\lambda)",
        xticks=(0:100:400, [L"%$s" for s in 0:100:400]),
        yticks=(-3:0, [L"%$s" for s in -3:0]),
        title=L"(b): Profiles for each $u_0$",
        titlealign=:left,
        width=600,
        height=400)
    for (i, result) in enumerate(profile_results)
        #prof = result.Cprof[:Dλ]
        #θ = get_parameter_values(prof)
        #ℓ = get_profile_values(prof)
        #mle_val = prof.parent.likelihood_solution[:Dλ]
        θ = result.θ
        ℓ = result.ℓ
        mle_val = result.mle_val
        lines!(ax2, θ, ℓ, color=colors[i], linewidth=3)
        linesegments!(ax2, [(mle_val, threshold), (mle_val, 0.0)], color=colors[i], linewidth=3, linestyle=:dash)
    end
    hlines!(ax2, [threshold], color=:red)
    ylims!(ax2, threshold - 1, 0.1)
    xlims!(ax2, 0, 400)

    ax3 = Axis(fig1[1, 2],
        xlabel=L"$t$ [day]",
        ylabel=L"\mu_{\mathrm{c}}(t)",
        title=L"(c):$ $ Square coverage",
        titlealign=:left,
        width=600,
        height=400,
        xticks=(0:10:40, [L"%$s" for s in 0:10:40]),
        yticks=(0:0.25:1, [L"%$s" for s in 0:0.25:1]))
    for (i, result) in enumerate(all_results)
        t = result.saveat_times_fine
        nt = length(t)
        idx_range = 1:nt
        intervals = result.intervals[idx_range]
        ℓ = first.(intervals)
        u = last.(intervals)
        mle_curve = result.q_mle
        m = mle_curve[idx_range]
        lines!(ax3, t, ℓ, linewidth=linewidth, color=colors[i])
        lines!(ax3, t, u, linewidth=linewidth, color=colors[i])
        lines!(ax3, t, m, linewidth=linewidth, linestyle=:dash, color=colors[i])
    end
    scatter!(ax3, reduce(vcat, days), reduce(vcat, areas), color=data_markercolor, markersize=markersize)
    xlims!(ax3, 0, 42)
    ylims!(ax3, 0, 1)

    ax4 = Axis(fig1[1, 3],
        xlabel=L"$t$ [day]",
        ylabel=L"$\mu_{\mathrm{p}}(t)",
        title=L"(e):$ $ Square perimeter",
        titlealign=:left,
        width=600,
        height=400,
        xticks=(0:10:40, [L"%$s" for s in 0:10:40]),
        yticks=(0:0.25:1, [L"%$s" for s in 0:0.25:1]))
    for (i, result) in enumerate(all_results)
        t = result.saveat_times_fine
        nt = length(t)
        idx_range = (nt+1):(2nt)
        intervals = result.intervals[idx_range]
        ℓ = first.(intervals)
        u = last.(intervals)
        mle_curve = result.q_mle
        m = mle_curve[idx_range]
        lines!(ax4, t, ℓ, linewidth=linewidth, color=colors[i])
        lines!(ax4, t, u, linewidth=linewidth, color=colors[i])
        lines!(ax4, t, m, linewidth=linewidth, linestyle=:dash, color=colors[i])
    end
    scatter!(ax4, reduce(vcat, days), reduce(vcat, perimeters), color=data_markercolor, markersize=markersize)
    xlims!(ax4, 0, 42)
    ylims!(ax4, 0, 1)

    ax5 = Axis(fig1[1, 4],
        xlabel=L"$t_b$ [day]",
        ylabel=L"p(t_b)",
        title=L"(g): Square $p(t_b)$",
        titlealign=:left,
        width=600,
        height=400,
        xticks=square_ptb_xticks,
        yticks=(0:0.1:0.2, [L"%$s" for s in 0:0.1:0.2]))
    for (i, result) in enumerate(all_results)
        q_vals = result.q_vals
        tb_vals = q_vals[end, :]
        t = result.saveat_times_fine
        nt = length(t)
        tb_idx = 3nt + 1
        mle_curve = result.q_mle
        tb_mle = mle_curve[tb_idx]
        density!(ax5, tb_vals, strokecolor=colors[i], color=(colors[i], 0.0), strokewidth=strokewidth, strokearound=strokearound)
        vlines!(ax5, [tb_mle], linestyle=:dash, linewidth=linewidth, color=colors[i])
    end

    ax6 = Axis(fig1[2, 2],
        xlabel=L"$t$ [day]",
        ylabel=L"$\mu_{\mathrm{c}}(t)",
        title=L"(d):$ $ Wave coverage",
        titlealign=:left,
        width=600,
        height=400,
        xticks=(0:10:40, [L"%$s" for s in 0:10:40]),
        yticks=(0:0.25:1, [L"%$s" for s in 0:0.25:1]))
    for (i, result) in enumerate(all_results)
        t = result.wave_saveat_times_fine
        nt = length(t)
        idx_range = 1:nt
        intervals = result.wave_intervals[idx_range]
        ℓ = first.(intervals)
        u = last.(intervals)
        mle_curve = result.wave_q_mle
        m = mle_curve[idx_range]
        lines!(ax6, t, ℓ, linewidth=linewidth, color=colors[i])
        lines!(ax6, t, u, linewidth=linewidth, color=colors[i])
        lines!(ax6, t, m, linewidth=linewidth, linestyle=:dash, color=colors[i])
    end
    scatter!(ax6, wave_days, wave_areas, color=data_markercolor, markersize=markersize)
    xlims!(ax6, 0, 42)
    ylims!(ax6, 0, 1)

    ax7 = Axis(fig1[2, 3],
        xlabel=L"$t$ [day]",
        ylabel=L"$\mu_{\mathrm{p}}(t)",
        title=L"(f):$ $ Wave perimeter",
        titlealign=:left,
        width=600,
        height=400,
        xticks=(0:10:40, [L"%$s" for s in 0:10:40]),
        yticks=(0:0.25:1, [L"%$s" for s in 0:0.25:1]))
    for (i, result) in enumerate(all_results)
        t = result.wave_saveat_times_fine
        nt = length(t)
        idx_range = (nt+1):(2nt)
        intervals = result.wave_intervals[idx_range]
        ℓ = first.(intervals)
        u = last.(intervals)
        mle_curve = result.wave_q_mle
        m = mle_curve[idx_range]
        lines!(ax7, t, ℓ, linewidth=linewidth, color=colors[i])
        lines!(ax7, t, u, linewidth=linewidth, color=colors[i])
        lines!(ax7, t, m, linewidth=linewidth, linestyle=:dash, color=colors[i])
    end
    scatter!(ax7, wave_days, wave_perimeters, color=data_markercolor, markersize=markersize)
    xlims!(ax7, 0, 42)
    ylims!(ax7, 0, 1)

    ax8 = Axis(fig1[2, 4],
        xlabel=L"$t_b$ [day]",
        ylabel=L"p(t_b)",
        title=L"(h): Wave $p(t_b)$",
        titlealign=:left,
        width=600,
        height=400,
        xticks=wave_ptb_xticks,
        yticks=(0:0.1:0.3, [L"%$s" for s in 0:0.1:0.3]))
    for (i, result) in enumerate(all_results)
        q_vals = result.wave_q_vals
        tb_vals = q_vals[end, :]
        t = result.wave_saveat_times_fine
        nt = length(t)
        tb_idx = 3nt + 1
        mle_curve = result.wave_q_mle
        tb_mle = mle_curve[tb_idx]
        density!(ax8, tb_vals, strokecolor=colors[i], color=(colors[i], 0.0), strokewidth=strokewidth, strokearound=strokearound)
        vlines!(ax8, [tb_mle], linestyle=:dash, linewidth=linewidth, color=colors[i])
    end

    Legend(fig1[1:2, 5], [LineElement(color=colors[i], linestyle=:solid, linewidth=9) for i in eachindex(u0)],
        [L"%$s" for s in u0],
        L"u_0", labelsize=64, titlesize=64, orientation=:vertical, titleposition=:top, halign=:left)
    resize_to_layout!(fig1)
    fig1
end

function plot_wave_likelihood_analysis(all_results, u0, pore_data, save_options, profile_results;
    Dλ_lim=(0.0, 600.0),
    square_ptb_xticks=(15:4:28, [L"%$s" for s in 15:4:28]),
    wave_ptb_xticks=(12:4:34, [L"%$s" for s in 12:4:34])
)
    colors = [:black, :blue, :magenta]
    levels = [ProfileLikelihood.get_chisq_threshold(0.95, 2)]
    threshold = ProfileLikelihood.get_chisq_threshold(0.95, 1)
    linewidth = 4
    markersize = 27
    data_markercolor = :blue
    strokearound = false
    strokewidth = 4
    areas, perimeters, radii, σ = get_pore_data(pore_data.data, save_options.saveat_times; combine=true)
    wave_areas, wave_perimeters, wave_radii, wave_σ = get_pore_data(pore_data.data, save_options.saveat_times[2], "SWV500"; combine=true)
    days = [repeat([save_options.saveat_times[i]], length(areas[i])) for i in eachindex(save_options.saveat_times)]
    wave_days = repeat([save_options.saveat_times[2]], length(wave_areas))

    fig1 = Figure(fontsize=63)
    ax1 = Axis(fig1[1, 1],
        xlabel=L"$D\lambda$ [$\mu$m$^2$/day$^2$]",
        ylabel=L"$\lambda$ [1/day]",
        title=L"(a):$ $ Likelihood CRs",
        titlealign=:left,
        xticks=(0:200:600, [L"%$s" for s in 0:200:600]),
        yticks=(0:2.5:5, [L"%$s" for s in 0:2.5:5.0]),
        width=600,
        height=400)
    for (i, result) in enumerate(all_results)
        grid_1 = result.param_ranges[1]
        grid_2 = result.param_ranges[2]
        normalised_loglik_vals_product = result.loglik_vals_product
        contour!(ax1, grid_1, grid_2, normalised_loglik_vals_product, levels=levels, color=colors[i], linewidth=linewidth)
        D̂ = profile_results[i].mle_val
        vlines!(ax1, [D̂], linestyle=:dash, linewidth=linewidth, color=colors[i])
    end
    xlims!(ax1, Dλ_lim...)

    ax2 = Axis(fig1[2, 1],
        xlabel=L"$D\lambda$ [$\mu$m$^2$/day$^2$]",
        ylabel=L"\ell_p(D\lambda)",
        xticks=(0:200:600, [L"%$s" for s in 0:200:600]),
        yticks=(-3:0, [L"%$s" for s in -3:0]),
        title=L"(b): Profiles for each $u_0$",
        titlealign=:left,
        width=600,
        height=400)
    for (i, result) in enumerate(profile_results)
        #prof = result.Cprof[:Dλ]
        #θ = get_parameter_values(prof)
        #ℓ = get_profile_values(prof)
        #mle_val = prof.parent.likelihood_solution[:Dλ]
        θ = result.θ
        ℓ = result.ℓ
        mle_val = result.mle_val
        lines!(ax2, θ, ℓ, color=colors[i], linewidth=3)
        linesegments!(ax2, [(mle_val, threshold), (mle_val, 0.0)], color=colors[i], linewidth=3, linestyle=:dash)
    end
    hlines!(ax2, [threshold], color=:red)
    ylims!(ax2, threshold - 1, 0.1)
    xlims!(ax2, Dλ_lim...)

    ax3 = Axis(fig1[1, 2],
        xlabel=L"$t$ [day]",
        ylabel=L"\mu_{\mathrm{c}}(t)",
        title=L"(c):$ $ Wave coverage",
        titlealign=:left,
        width=600,
        height=400,
        xticks=(0:10:40, [L"%$s" for s in 0:10:40]),
        yticks=(0:0.25:1, [L"%$s" for s in 0:0.25:1]))
    for (i, result) in enumerate(all_results)
        t = result.saveat_times_fine
        nt = length(t)
        idx_range = 1:nt
        intervals = result.intervals[idx_range]
        ℓ = first.(intervals)
        u = last.(intervals)
        mle_curve = result.q_mle
        m = mle_curve[idx_range]
        lines!(ax3, t, ℓ, linewidth=linewidth, color=colors[i])
        lines!(ax3, t, u, linewidth=linewidth, color=colors[i])
        lines!(ax3, t, m, linewidth=linewidth, linestyle=:dash, color=colors[i])
    end
    scatter!(ax3, reduce(vcat, wave_days), reduce(vcat, wave_areas), color=data_markercolor, markersize=markersize)
    xlims!(ax3, 0, 42)
    ylims!(ax3, 0, 1)

    ax4 = Axis(fig1[1, 3],
        xlabel=L"$t$ [day]",
        ylabel=L"$\mu_{\mathrm{p}}(t)",
        title=L"(e):$ $ Wave perimeter",
        titlealign=:left,
        width=600,
        height=400,
        xticks=(0:10:40, [L"%$s" for s in 0:10:40]),
        yticks=(0:0.25:1, [L"%$s" for s in 0:0.25:1]))
    for (i, result) in enumerate(all_results)
        t = result.saveat_times_fine
        nt = length(t)
        idx_range = (nt+1):(2nt)
        intervals = result.intervals[idx_range]
        ℓ = first.(intervals)
        u = last.(intervals)
        mle_curve = result.q_mle
        m = mle_curve[idx_range]
        lines!(ax4, t, ℓ, linewidth=linewidth, color=colors[i])
        lines!(ax4, t, u, linewidth=linewidth, color=colors[i])
        lines!(ax4, t, m, linewidth=linewidth, linestyle=:dash, color=colors[i])
    end
    scatter!(ax4, reduce(vcat, wave_days), reduce(vcat, wave_perimeters), color=data_markercolor, markersize=markersize)
    xlims!(ax4, 0, 42)
    ylims!(ax4, 0, 1)

    ax5 = Axis(fig1[1, 4],
        xlabel=L"$t_b$ [day]",
        ylabel=L"p(t_b)",
        title=L"(g): Wave $p(t_b)$",
        titlealign=:left,
        width=600,
        height=400,
        xticks=wave_ptb_xticks,
        yticks=(0:0.1:0.2, [L"%$s" for s in 0:0.1:0.2]))
    for (i, result) in enumerate(all_results)
        q_vals = result.q_vals
        tb_vals = q_vals[end, :]
        t = result.saveat_times_fine
        nt = length(t)
        tb_idx = 3nt + 1
        mle_curve = result.q_mle
        tb_mle = mle_curve[tb_idx]
        density!(ax5, tb_vals, strokecolor=colors[i], color=(colors[i], 0.0), strokewidth=strokewidth, strokearound=strokearound)
        vlines!(ax5, [tb_mle], linestyle=:dash, linewidth=linewidth, color=colors[i])
    end

    ax6 = Axis(fig1[2, 2],
        xlabel=L"$t$ [day]",
        ylabel=L"$\mu_{\mathrm{c}}(t)",
        title=L"(d):$ $ Square coverage",
        titlealign=:left,
        width=600,
        height=400,
        xticks=(0:10:40, [L"%$s" for s in 0:10:40]),
        yticks=(0:0.25:1, [L"%$s" for s in 0:0.25:1]))
    for (i, result) in enumerate(all_results)
        t = result.square_saveat_times_fine
        nt = length(t)
        idx_range = 1:nt
        intervals = result.square_intervals[idx_range]
        ℓ = first.(intervals)
        u = last.(intervals)
        mle_curve = result.square_q_mle
        m = mle_curve[idx_range]
        lines!(ax6, t, ℓ, linewidth=linewidth, color=colors[i])
        lines!(ax6, t, u, linewidth=linewidth, color=colors[i])
        lines!(ax6, t, m, linewidth=linewidth, linestyle=:dash, color=colors[i])
    end
    scatter!(ax6, reduce(vcat, days), reduce(vcat, areas), color=data_markercolor, markersize=markersize)
    xlims!(ax6, 0, 42)
    ylims!(ax6, 0, 1)

    ax7 = Axis(fig1[2, 3],
        xlabel=L"$t$ [day]",
        ylabel=L"$\mu_{\mathrm{p}}(t)",
        title=L"(f):$ $ Square perimeter",
        titlealign=:left,
        width=600,
        height=400,
        xticks=(0:10:40, [L"%$s" for s in 0:10:40]),
        yticks=(0:0.25:1, [L"%$s" for s in 0:0.25:1]))
    for (i, result) in enumerate(all_results)
        t = result.square_saveat_times_fine
        nt = length(t)
        idx_range = (nt+1):(2nt)
        intervals = result.square_intervals[idx_range]
        ℓ = first.(intervals)
        u = last.(intervals)
        mle_curve = result.square_q_mle
        m = mle_curve[idx_range]
        lines!(ax7, t, ℓ, linewidth=linewidth, color=colors[i])
        lines!(ax7, t, u, linewidth=linewidth, color=colors[i])
        lines!(ax7, t, m, linewidth=linewidth, linestyle=:dash, color=colors[i])
    end
    scatter!(ax7, reduce(vcat, days), reduce(vcat, perimeters), color=data_markercolor, markersize=markersize)
    xlims!(ax7, 0, 42)
    ylims!(ax7, 0, 1)

    ax8 = Axis(fig1[2, 4],
        xlabel=L"$t_b$ [day]",
        ylabel=L"p(t_b)",
        title=L"(h): Square $p(t_b)$",
        titlealign=:left,
        width=600,
        height=400,
        xticks=square_ptb_xticks,
        yticks=(0:0.1:0.3, [L"%$s" for s in 0:0.1:0.3]))
    for (i, result) in enumerate(all_results)
        q_vals = result.square_q_vals
        tb_vals = q_vals[end, :]
        t = result.square_saveat_times_fine
        nt = length(t)
        tb_idx = 3nt + 1
        mle_curve = result.square_q_mle
        tb_mle = mle_curve[tb_idx]
        density!(ax8, tb_vals, strokecolor=colors[i], color=(colors[i], 0.0), strokewidth=strokewidth, strokearound=strokearound)
        vlines!(ax8, [tb_mle], linestyle=:dash, linewidth=linewidth, color=colors[i])
    end

    Legend(fig1[1:2, 5], [LineElement(color=colors[i], linestyle=:solid, linewidth=9) for i in eachindex(u0)],
        [L"%$s" for s in u0],
        L"u_0", labelsize=64, titlesize=64, orientation=:vertical, titleposition=:top, halign=:left)
    resize_to_layout!(fig1)
    fig1
end

function plot_tissue_growth_predictions(result, pore_data, save_options)
    # square 
    mesh = pore_data.mesh
    snapshot_save_times = save_options.snapshot_save_times
    centroid = pore_data.centroid
    _D = result.C_D
    _λ = result.C_λ
    _u₀ = result.C_u₀
    lad1, lad2, sols = get_all_leading_edges(_D, _λ, _u₀, mesh, snapshot_save_times, centroid)
    pts = mesh.mesh_information.triangulation.points
    tri = mesh.mesh_information.triangulation
    bn = tri.boundary_nodes
    pts_bn = pts[:, bn]
    b1 = pts_bn[1, :]
    b2 = pts_bn[2, :]
    fig1 = plot_tissue_growth_predictions(lad1, lad2, sols, b1, b2, pts, tri, snapshot_save_times)

    # wave
    mesh = pore_data.wave_mesh
    snapshot_save_times = save_options.snapshot_save_times
    centroid = pore_data.wave_centroid
    lad1, lad2, sols = get_all_leading_edges(_D, _λ, _u₀, mesh, snapshot_save_times, centroid; use_concave_hull=true, rotate=true)
    pts = mesh.mesh_information.triangulation.points
    tri = mesh.mesh_information.triangulation
    bn = tri.boundary_nodes
    pts_bn = pts[:, bn]
    b1 = pts_bn[2, :]
    b2 = pts_bn[1, :]
    pts = [pts[2, :]'; pts[1, :]']
    fig2 = plot_tissue_growth_predictions(lad1, lad2, sols, b1, b2, pts, tri, snapshot_save_times)
    return fig1, fig2
end
function plot_tissue_growth_predictions(lad1, lad2, sols, b1, b2, pts, tri, snapshot_save_times)
    fig1 = Figure(fontsize=63, resolution=(3255, 2330))
    alph_index = 0
    for i in 1:3
        for (j, τ) in pairs(snapshot_save_times)
            if j > 1
                alph_index += 1
                if i == 1
                    ax = Axis(
                        fig1[i, j-1],
                        xlabel=L"$x$ [$\mu$m$]",
                        ylabel=L"$y$ [$\mu$m$]",
                        title=L"$ $Day %$(τ)",
                        titlealign=:left,
                        aspect=1,
                        width=600,
                        height=600,
                        xticks=(0:200:600, [L"%$s" for s in 0:200:600]),
                        yticks=(0:200:600, [L"%$s" for s in 0:200:600]),
                        titlesize=88
                    )
                else
                    ax = Axis(
                        fig1[i, j-1],
                        xlabel=L"$x$ [$\mu$m$]",
                        ylabel=L"$y$ [$\mu$m$]",
                        aspect=1,
                        width=600,
                        height=600,
                        xticks=(0:200:600, [L"%$s" for s in 0:200:600]),
                        yticks=(0:200:600, [L"%$s" for s in 0:200:600])
                    )
                end
                tricontourf!(ax, pts[1, :], pts[2, :], sols[i][j], levels=range(0.48, 1.1, length=25), extendlow=:black, extendhigh=:green, colormap=:algae, triangulation=[T[j] for T in each_triangle(tri), j in 1:3]')
                if !isnothing(lad1[i][j])
                    lines!(ax, lad1[i][j], lad2[i][j], color=:red, linewidth=9)
                end
                lines!(ax, b1, b2, linewidth=3, color=:black)
                if j > 2
                    hideydecorations!(ax; grid=false)
                end
                if i < 3
                    hidexdecorations!(ax; grid=false)
                end
                tightlimits!(ax)
            end
        end
    end
    Label(fig1[1, 0], L"\hat{\mathbf{\theta}}_L", fontsize=88)
    Label(fig1[2, 0], L"\hat{\mathbf{\theta}}", fontsize=88)
    Label(fig1[3, 0], L"\hat{\mathbf{\theta}}_U", fontsize=88)
    resize_to_layout!(fig1)
    return fig1
end