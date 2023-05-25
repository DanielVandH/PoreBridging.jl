"""
    pore_prediction_function(θ::AbstractVector{TT}, p) where {TT}

Returns the prediction function for the given parameters `θ = (D, λ, u₀)` and parameters `p`. See 
[`compute_prediction_intervals`](@ref).
"""
function pore_prediction_function(θ::AbstractVector{TT}, p) where {TT}
    D, λ, u₀ = θ
    (; mesh, ext_initial_time, ext_final_time, saveat_times_fine,
        pore_area, pore_perimeter, pore_smallest_rad, centroid, use_concave_hull) = p
    BC = get_boundary_conditions(mesh, λ)
    initial_time = ext_initial_time
    final_time = ext_final_time
    prob, jac_prototype = get_pde(mesh, BC, D, λ, u₀, initial_time, final_time)
    cb = ContinuousCallback(
        (u, t, integrator) -> begin
            T = FVM.get_elements(prob)
            pts = FVM.get_points(prob)
            A, _, _ = compute_leading_edge_statistics(u, T, pts, 1 / 2, centroid)
            A = 1.0 - A / pore_area
            return A - 1e-9
        end,
        terminate!
    )
    sol = solve(prob,
        TRBDF2(linsolve=KLUFactorization(; reuse_symbolic=false));
        jac_prototype=float.(jac_prototype),
        saveat=saveat_times_fine,
        callback=cb,
        parallel=true)
    n = length(saveat_times_fine)
    q = zeros(3n + 1)
    area_cache = zeros(n)
    perimeter_cache = zeros(n)
    radius_cache = zeros(n)
    closing_time = sol.t[end]
    for i in eachindex(saveat_times_fine)
        if i ≤ length(sol) && (sol.t[i] ∈ saveat_times_fine)
            if !use_concave_hull
                A, P, R = compute_leading_edge_statistics(
                    sol.u[i],
                    FVM.get_elements(prob),
                    FVM.get_points(prob),
                    1 / 2,
                    centroid
                )
                area_cache[i] = 1.0 - A / pore_area
                perimeter_cache[i] = iszero(P) ? 1.0 : P / pore_perimeter
                radius_cache[i] = isinf(R) ? 1.0 : R / pore_smallest_rad
            else
                _u = sol.u[i]
                lad = get_leading_edge(_u, prob, 1 / 2, centroid)
                if !isempty(lad)
                    append!(lad, lad[:, begin])
                    ch = concave_hull([[x, y] for (x, y) in eachcol(lad)])
                    xy_bn = ch.vertices
                    push!(xy_bn, xy_bn[begin])
                    x_bn = first.(xy_bn)
                    y_bn = last.(xy_bn)
                    area_cache[i] = abs(PolygonOps.area(xy_bn) / pore_area)
                    if iszero(area_cache[i]) && sol.t[i] < 14.0
                        area_cache[i] = 1.0
                    end
                    perimeter_cache[i] = sum(sqrt((x_bn[i+1] - x_bn[i])^2 + (y_bn[i+1] - y_bn[i])^2) for i in 1:(length(x_bn)-1)) / pore_perimeter
                    if iszero(perimeter_cache[i]) && sol.t[i] < 14.0
                        perimeter_cache[i] = 1.0
                    end
                    xy_mat = [x_bn'; y_bn']
                    radius_cache[i] = get_smallest_distance(xy_mat, centroid) / pore_smallest_rad
                    if iszero(radius_cache[i]) || isinf(radius_cache[i]) && sol.t[i] < 14.0
                        radius_cache[i] = 1.0
                    end
                end
            end
        else
            area_cache[i] = 0.0
            perimeter_cache[i] = 0.0
            radius_cache[i] = 0.0
        end
    end
    q[1:n] .= area_cache
    q[(n+1):2n] .= perimeter_cache
    q[(2n+1):3n] .= radius_cache
    q[3n+1] = closing_time
    @show closing_time
    return q
end

"""
    construct_prediction_function(likprob)

Returns a mapping frmo the parameters in `likprob` to the form `(D, λ, u₀)` used in 
the prediction function.
"""
function construct_prediction_function(likprob)
    fnc = (θ, p) -> pore_prediction_function(likprob.log_likelihood_function.θ_conv.contents(θ), p)
    return fnc
end

"""
    compute_prediction_intervals(prof,
        nt=100,
        initial_time=prof.likelihood_problem.data.initial_time,
        final_time=5prof.likelihood_problem.data.final_time;
        resolution=50,
        parallel=true,
        mesh=prof.likelihood_problem.data.mesh,
        ext_initial_time=initial_time,
        ext_final_time=final_time,
        pore_area=prof.likelihood_problem.data.pore_area,
        pore_perimeter=prof.likelihood_problem.data.pore_perimeter,
        pore_smallest_rad=prof.likelihood_problem.data.pore_smallest_rad,
        centroid=prof.likelihood_problem.data.centroid,
        use_concave_hull=false)

Compute the prediction intervals from the given `ProfileLikelihoodSolution`, `prof`.

# Arguments
- `prof`: The `ProfileLikelihoodSolution`.
- `nt=100`: The number of time points to compute predictions at.
- `initial_time=prof.likelihood_problem.data.initial_time`: The initial time for the model.
- `final_time=5prof.likelihood_problem.data.final_time`: The final time for the model. This should be sufficiently large so that the bridging time can be found in all cases.
- `resolution=50`: The number of points to use for the parameters.
- `parallel=true`: Whether to use multithreading.
- `mesh=prof.likelihood_problem.data.mesh`: The mesh.
- `ext_initial_time=initial_time`: The initial time for the prediction function.
- `ext_final_time=final_time`: The final time for the prediction function.
- `pore_area=prof.likelihood_problem.data.pore_area`: The pore area.
- `pore_perimeter=prof.likelihood_problem.data.pore_perimeter`: The pore perimeter.
- `pore_smallest_rad=prof.likelihood_problem.data.pore_smallest_rad`: The smallest distance from the `centroid` to the pore boundary.
- `centroid=prof.likelihood_problem.data.centroid`: The centroid of the pore.
- `use_concave_hull=false`: Whether to use the concave hull for the prediction function.

# Outputs
- `individual_intervals`: The individual intervals for each parameter for the prediction intervals.
- `union_intervals`: The unions of the individual intervals for each parameter for the prediction intervals.
- `q_vals`: The values of `q` for the prediction intervals.
- `param_ranges`: The parameter ranges for the prediction intervals.
- `q_mle`: The value of the prediction function at the MLE.
- `saveat_times_fine`: The times to save the numerical solutions at for the objective function when considering the prediction intervals.
"""
function compute_prediction_intervals(prof,
    nt=100,
    initial_time=prof.likelihood_problem.data.initial_time,
    final_time=5prof.likelihood_problem.data.final_time;
    resolution=50,
    parallel=true,
    mesh=prof.likelihood_problem.data.mesh,
    ext_initial_time=initial_time,
    ext_final_time=final_time,
    pore_area=prof.likelihood_problem.data.pore_area,
    pore_perimeter=prof.likelihood_problem.data.pore_perimeter,
    pore_smallest_rad=prof.likelihood_problem.data.pore_smallest_rad,
    centroid=prof.likelihood_problem.data.centroid,
    use_concave_hull=false)
    likprob = prof.likelihood_problem
    mle_sol = prof.likelihood_solution
    saveat_times_fine = LinRange(ext_initial_time, ext_final_time, nt)
    p = (mesh=mesh,
        ext_initial_time=ext_initial_time,
        ext_final_time=ext_final_time,
        saveat_times_fine=saveat_times_fine,
        pore_area=pore_area,
        pore_perimeter=pore_perimeter,
        pore_smallest_rad=pore_smallest_rad,
        centroid=centroid,
        use_concave_hull=use_concave_hull)
    fnc = construct_prediction_function(likprob)
    q_prototype = zeros(3nt + 1)
    individual_intervals, union_intervals, q_vals, param_ranges =
        get_prediction_intervals(
            fnc,
            prof,
            p;
            parallel,
            q_prototype,
            resolution
        )
    q_mle = fnc(get_mle(mle_sol), p)
    return individual_intervals, union_intervals, q_vals, param_ranges, q_mle, saveat_times_fine
end

"""
    plot_prediction_functions(saveat_times_fine,
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
        fix_u₀
    )

Plot the prediction intervals. See also [`compute_prediction_intervals`](@ref).

# Arguments
- `saveat_times_fine`: The times to save the numerical solutions at for the objective function when considering the prediction intervals.
- `individual_intervals`: The individual intervals for each parameter for the prediction intervals.
- `union_intervals`: The unions of the individual intervals for each parameter for the prediction intervals.
- `q_mle`: The value of the prediction function at the MLE.
- `areas`: The void coverages.
- `perimeters`: The normalised void perimeters.
- `radii`: The normalised void radii.
- `saveat_times`: The times to save the numerical solutions at for the objective function.
- `product`: Whether to consider `Dλ` or `(D, λ)`.
- `fix_D`: Whether to fix `D`.
- `fix_λ`: Whether to fix `λ`.
- `fix_u₀`: Whether to fix `u₀`.

# Outputs
- `fig`: The figure of the prediction intervals.
- `indiv_idx`: The indices of the individual intervals for each parameter for the prediction intervals.
"""
function plot_prediction_functions(saveat_times_fine,
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
    fix_u₀
)
    if !product
        if !fix_D && !fix_λ && !fix_u₀
            titles = (
                alph -> L"(%$(alph)): Profile-wise PI for $D$",
                alph -> L"(%$(alph)): Profile-wise PI for $\lambda$",
                alph -> L"(%$(alph)): Profile-wise PI for $u_0$",
                alph -> L"(%$(alph)):$ $ Union of all intervals"
            )
            alph = ['a' 'b' 'c' 'd'
                'e' 'f' 'g' 'h'
                'i' 'j' 'k' 'ℓ']
            lower_union = first.(union_intervals)
            upper_union = last.(union_intervals)

            indiv_idx = [1, 2, 3]
        elseif !fix_D && !fix_λ && fix_u₀
            titles = (
                alph -> L"(%$(alph)): Profile-wise PI for $D$",
                alph -> L"(%$(alph)): Profile-wise PI for $\lambda$",
                alph -> L"(%$(alph)):$ $ Union of all intervals"
            )
            alph = ['a' 'b' 'c'
                'd' 'e' 'f'
                'g' 'h' 'i']
            lower_union = first.(union_intervals)
            upper_union = last.(union_intervals)
            indiv_idx = [1, 2]
        elseif !fix_D && fix_λ && !fix_u₀
            titles = (
                alph -> L"(%$(alph)): Profile-wise PI for $D$",
                alph -> L"(%$(alph)): Profile-wise PI for $u_0$",
                alph -> L"(%$(alph)):$ $ Union of all intervals"
            )
            alph = ['a' 'b' 'c'
                'd' 'e' 'f'
                'g' 'h' 'i']
            lower_union = first.(union_intervals)
            upper_union = last.(union_intervals)
            indiv_idx = [1, 2]
        elseif fix_D && !fix_λ && !fix_u₀
            titles = (
                alph -> L"(%$(alph)): Profile-wise PI for $\lambda$",
                alph -> L"(%$(alph)): Profile-wise PI for $u_0$",
                alph -> L"(%$(alph)):$ $ Union of all intervals"
            )
            alph = ['a' 'b' 'c'
                'd' 'e' 'f'
                'g' 'h' 'i']
            lower_union = first.(union_intervals)
            upper_union = last.(union_intervals)
            indiv_idx = [1, 2]
        elseif !fix_D && fix_λ && fix_u₀
            titles = (
                alph -> L"(%$(alph)): Profile-wise PI for $D$",
                nothing
            )
            alph = ['a'
                'b'
                'c']
            lower_union = first.(union_intervals)
            upper_union = last.(union_intervals)
            indiv_idx = [1]
        elseif fix_D && !fix_λ && fix_u₀
            titles = (
                alph -> L"(%$(alph)): Profile-wise PI for $\lambda$",
                nothing
            )
            alph = ['a'
                'b'
                'c']
            lower_union = first.(union_intervals)
            upper_union = last.(union_intervals)
            indiv_idx = [1]
        elseif fix_D && fix_λ && !fix_u₀
            titles = (
                alph -> L"(%$(alph)): Profile-wise PI for $u_0$",
                nothing
            )
            alph = ['a'
                'b'
                'c']
            lower_union = first.(union_intervals)
            upper_union = last.(union_intervals)
            indiv_idx = [1]
        elseif fix_D && fix_λ && fix_u₀
            throw("Cannot fix all three parameters.")
        end
    elseif product
        if !fix_D && !fix_λ && !fix_u₀
            titles = (
                alph -> L"(%$(alph)): Profile-wise PI for $D\lambda$",
                alph -> L"(%$(alph)): Profile-wise PI for $u_0$",
                alph -> L"(%$(alph)):$ $ Union of all intervals"
            )
            alph = ['a' 'b' 'c'
                'e' 'f' 'g'
                'i' 'j' 'k']
            lower_union = [min(a, b) for (a, b) in zip(first.(individual_intervals[1]), first.(individual_intervals[3]))]
            upper_union = [max(a, b) for (a, b) in zip(last.(individual_intervals[1]), last.(individual_intervals[3]))]
            indiv_idx = [1, 3]
        elseif !fix_D && !fix_λ && fix_u₀
            titles = (
                alph -> L"(%$(alph)): Profile-wise PI for $D\lambda$",
                nothing
            )
            alph = ['a'
                'b'
                'c']
            lower_union = first.(individual_intervals[1])
            upper_union = last.(individual_intervals[1])
            indiv_idx = [1]
        elseif !fix_D && fix_λ && !fix_u₀
            throw("Cannot fix λ while optimising Dλ.")
        elseif fix_D && !fix_λ && !fix_u₀
            throw("Cannot fix D while optimising Dλ")
        elseif !fix_D && fix_λ && fix_u₀
            throw("Cannot fix λ while optimsing Dλ.")
        elseif fix_D && !fix_λ && fix_u₀
            throw("Cannot fix D while optimising Dλ")
        elseif fix_D && fix_λ && !fix_u₀
            throw("Cannot fix D and λ while optimising Dλ")
        elseif fix_D && fix_λ && fix_u₀
            throw("Cannot fix all three parameters.")
        end
    end
    nt = length(saveat_times_fine)
    stat_indices = (1:nt, (nt+1):(2nt), (2nt+1):3nt)
    data_plot = [areas, perimeters, radii]
    days = [repeat([saveat_times[i]], length(areas[i])) for i in eachindex(saveat_times)]
    yaxes = (L"$ $Coverage", L"$ $Perimeter", L"$ $Radius")
    fig = Figure(fontsize=36, resolution=(1988, 1365))
    for i in axes(alph, 1)
        for j in axes(alph, 2)
            ax = Axis(
                fig[i, j],
                xlabel=L"$t$ (d)",
                ylabel=yaxes[i],
                title=titles[j](alph[i, j]),
                titlealign=:left,
                width=600,
                height=300,
                xticks=([7, 14, 21, 28, 35, 42], [L"%$s" for s in [7, 14, 21, 28, 35, 42]])
            )
            if j ≠ size(alph, 2) || isnothing(titles[end])
                band!(
                    ax,
                    saveat_times_fine,
                    first.(individual_intervals[indiv_idx[j]])[stat_indices[i]],
                    last.(individual_intervals[indiv_idx[j]])[stat_indices[i]],
                    color=(:blue, 0.35)
                )
            else
                band!(
                    ax,
                    saveat_times_fine,
                    lower_union[stat_indices[i]],
                    upper_union[stat_indices[i]],
                    color=(:blue, 0.35)
                )
            end
            lines!(
                ax,
                saveat_times_fine,
                q_mle[stat_indices[i]],
                linewidth=4,
                color=:red
            )
            j ≠ 1 && hideydecorations!(ax)
            ylims!(ax, 0, 1)
            xlims!(ax, 0, 42)
            scatter!(ax, reduce(vcat, days), reduce(vcat, data_plot[i]), color=:blue, markersize=27)
        end
    end
    return fig, indiv_idx
end

"""
    plot_closing_time_densities(q_vals, q_mle, indiv_idx)

Plot the closing time densities. See also [`compute_prediction_intervals`](@ref).

# Arguments
- `q_vals`: The values of `q` for the prediction intervals.
- `q_mle`: The value of the prediction function at the MLE.
- `indiv_idx`: The indices of the individual intervals for each parameter for the prediction intervals.

# Outputs
- `fig`: The figure of the closing time densities.
"""
function plot_closing_time_densities(q_vals, q_mle, indiv_idx)
    fig = Figure(fontsize=32)
    ax = Axis(
        fig[1, 1],
        xlabel=L"t_c",
        ylabel=L"$ $Probability density",
        width=600,
        height=300
    )
    density!(ax,
        reduce(vcat, [q_vals[i][end, :] for i in indiv_idx]),
        strokecolor=:blue,
        color=(:blue, 0.3),
        strokewidth=3,
        strokearound=false)
    vlines!(ax, [q_mle[end]], color=:red, linewidth=4)
    resize_to_layout!(fig)
    return fig
end