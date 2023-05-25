"""
    get_Dλu₀_for_snapshots(prof, mle_sol;
        product=false,
        fix_D=false,
        fix_λ=false,
        fix_u₀=false,
        fixed_D=nothing,
        fixed_λ=nothing,
        fixed_u₀=nothing)

Returns values of `(D, λ, u₀)` to use for the snapshots.

# Arguments
- `prof`: The `ProfileLikelihoodSolution`. See [`find_profiles`](@ref).
- `mle_sol`: The maximum likelihood solution. See [`find_mle`](@ref).

# Keyword Arguments
- `product=false`: Whether to consider `Dλ` or `(D, λ)`.
- `fix_D=false`: Whether to fix `D`.
- `fix_λ=false`: Whether to fix `λ`.
- `fix_u₀=false`: Whether to fix `u₀`.
- `fixed_D=nothing`: The fixed value for `D`, or `nothing` if `fix_D==false`.
- `fixed_λ=nothing`: The fixed value for `λ`, or `nothing` if `fix_λ==false`.
- `fixed_u₀=nothing`: The fixed value for `u₀`, or `nothing` if `fix_u₀==false`.

# Outputs 
- `_D`: The values of `D` to use for the snapshots, given by `(lower, mle, upper)`.
- `_λ`: The values of `λ` to use for the snapshots, given by `(lower, mle, upper)`.
- `_u₀`: The values of `u₀` to use for the snapshots, given by `(lower, mle, upper)`.

Here, `lower` and `upper` refer to the corresponding endpoints for the confidence intervals,
and `mle` is the estimate from `mle_sol`. Note that if `Dλ` is considered, then 
`Dλ/λ = D` is used. 
"""
function get_Dλu₀_for_snapshots(prof, mle_sol;
    product=false,
    fix_D=false,
    fix_λ=false,
    fix_u₀=false,
    fixed_D=nothing,
    fixed_λ=nothing,
    fixed_u₀=nothing)
    if !product
        if !fix_D && !fix_λ && !fix_u₀
            _D = [
                mle_sol[:D],
                get_confidence_intervals(prof, :D).lower,
                get_confidence_intervals(prof, :D).upper
            ]
            _λ = [
                mle_sol[:λ],
                get_confidence_intervals(prof, :λ).lower,
                get_confidence_intervals(prof, :λ).upper
            ]
            _u₀ = [
                mle_sol[:u₀],
                get_confidence_intervals(prof, :u₀).lower,
                get_confidence_intervals(prof, :u₀).upper
            ]
        elseif !fix_D && !fix_λ && fix_u₀
            _D = [
                mle_sol[:D],
                get_confidence_intervals(prof, :D).lower,
                get_confidence_intervals(prof, :D).upper
            ]
            _λ = [
                mle_sol[:λ],
                get_confidence_intervals(prof, :λ).lower,
                get_confidence_intervals(prof, :λ).upper
            ]
            _u₀ = [
                fixed_u₀,
                fixed_u₀,
                fixed_u₀
            ]
        elseif !fix_D && fix_λ && !fix_u₀
            _D = [
                mle_sol[:D],
                get_confidence_intervals(prof, :D).lower,
                get_confidence_intervals(prof, :D).upper
            ]
            _λ = [
                fixed_λ,
                fixed_λ,
                fixed_λ
            ]
            _u₀ = [
                mle_sol[:u₀],
                get_confidence_intervals(prof, :u₀).lower,
                get_confidence_intervals(prof, :u₀).upper
            ]
        elseif fix_D && !fix_λ && !fix_u₀
            _D = [
                fixed_D,
                fixed_D,
                fixed_D
            ]
            _λ = [
                mle_sol[:λ],
                get_confidence_intervals(prof, :λ).lower,
                get_confidence_intervals(prof, :λ).upper
            ]
            _u₀ = [
                mle_sol[:u₀],
                get_confidence_intervals(prof, :u₀).lower,
                get_confidence_intervals(prof, :u₀).upper
            ]
        elseif !fix_D && fix_λ && fix_u₀
            _D = [
                mle_sol[:D],
                get_confidence_intervals(prof, :D).lower,
                get_confidence_intervals(prof, :D).upper
            ]
            _λ = [
                fixed_λ,
                fixed_λ,
                fixed_λ
            ]
            _u₀ = [
                fixed_u₀,
                fixed_u₀,
                fixed_u₀
            ]
        elseif fix_D && !fix_λ && fix_u₀
            _D = [
                fixed_D,
                fixed_D,
                fixed_D
            ]
            _λ = [
                mle_sol[:λ],
                get_confidence_intervals(prof, :λ).lower,
                get_confidence_intervals(prof, :λ).upper
            ]
            _u₀ = [
                fixed_u₀,
                fixed_u₀,
                fixed_u₀
            ]
        elseif fix_D && fix_λ && !fix_u₀
            _D = [
                fixed_D,
                fixed_D,
                fixed_D
            ]
            _λ = [
                fixed_λ,
                fixed_λ,
                fixed_λ
            ]
            _u₀ = [
                mle_sol[:u₀],
                get_confidence_intervals(prof, :u₀).lower,
                get_confidence_intervals(prof, :u₀).upper
            ]
        elseif fix_D && fix_λ && fix_u₀
            throw("Cannot fix all three parameters.")
        end
    elseif product
        if !fix_D && !fix_λ && !fix_u₀
            _Dλ = [
                mle_sol[:Dλ],
                get_confidence_intervals(prof, :Dλ).lower,
                get_confidence_intervals(prof, :Dλ).upper
            ]
            _λ = [
                mle_sol[:λ],
                get_confidence_intervals(prof, :λ).lower,
                get_confidence_intervals(prof, :λ).upper
            ]
            _D = _Dλ ./ _λ
            _u₀ = [
                mle_sol[:u₀],
                get_confidence_intervals(prof, :u₀).lower,
                get_confidence_intervals(prof, :u₀).upper
            ]
        elseif !fix_D && !fix_λ && fix_u₀
            _Dλ = [
                mle_sol[:Dλ],
                get_confidence_intervals(prof, :Dλ).lower,
                get_confidence_intervals(prof, :Dλ).upper
            ]
            _λ = [
                mle_sol[:λ],
                get_confidence_intervals(prof, :λ).lower,
                get_confidence_intervals(prof, :λ).upper
            ]
            _D = _Dλ ./ _λ
            _u₀ = [
                fixed_u₀,
                fixed_u₀,
                fixed_u₀
            ]
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
    _D = _D[[2, 1, 3]]
    _λ = _λ[[2, 1, 3]]
    _u₀ = _u₀[[2, 1, 3]]
    return _D, _λ, _u₀
end

"""
    plot_snapshots(_D, _λ, _u₀, mesh, snapshot_save_times, centroid;
        use_concave_hull=false,
        rotate=false,
        shift_vals=false,
        show_text=false)

Plot the snapshots. See also [`get_Dλu₀_for_snapshots`](@ref).

# Arguments
- `_D`: The `D` values to use for the snapshots. 
- `_λ`: The `λ` values to use for the snapshots.
- `_u₀`: The `u₀` values to use for the snapshots.
- `mesh`: The mesh.
- `snapshot_save_times`: The times at which to save the snapshots.
- `centroid`: The centroid of the mesh. 

# Keyword Arguments
- `use_concave_hull=false`: Whether to use the concave hull of the leading edge for plotting.
- `rotate=false`: Whether to rotate the mesh by 90 degrees.
- `shift_vals=false`: Whether to shift the values to the top of the figure when showing the values with text.
- `show_text=false`: Whether to show the values with text.

# Outputs
- `fig`: The figure.
"""
function plot_snapshots(_D, _λ, _u₀, mesh, snapshot_save_times, centroid;
    use_concave_hull=false, rotate=false, shift_vals=false, show_text=false)
    fig = Figure(fontsize=43, resolution=(3255, 2330))
    initial_time = minimum(snapshot_save_times)
    final_time = maximum(snapshot_save_times)
    alph_index = 0
    for i in 1:3
        D = _D[i]
        λ = _λ[i]
        u₀ = _u₀[i]
        BC = get_boundary_conditions(mesh, λ)
        prob, jac_prototype = get_pde(mesh, BC, D, λ, u₀, initial_time, final_time)
        sol = solve(prob,
            TRBDF2(linsolve=KLUFactorization(; reuse_symbolic=true));
            jac_prototype=float.(jac_prototype),
            parallel=true,
            saveat=snapshot_save_times)
        for (j, τ) in pairs(snapshot_save_times)
            alph_index += 1
            ax = Axis(
                fig[i, j],
                xlabel=L"x",
                ylabel=L"y",
                title=L"(%$(join('a':'z')[alph_index])): $t = %$(τ)$ d",
                titlealign=:left,
                aspect=1,
                width=600,
                height=600,
                xticks=(0:200:600, [L"%$s" for s in 0:200:600]),
                yticks=(0:200:600, [L"%$s" for s in 0:200:600])
            )
            pts = FVM.get_points(prob)
            tri = prob.mesh.mesh_information.triangulation
            if !rotate
                mesh!(ax, pts, [T[j] for T in each_triangle(tri), j in 1:3], color=sol.u[j], colorrange=(0.0, 1.0), colormap=:binary)
                lines!(ax, pts[:, get_boundary_nodes(tri)], linewidth=3, color=:black)
            else
                p1 = pts[1, :]
                p2 = pts[2, :]
                mesh!(ax, [p2'; p1'], [T[j] for T in each_triangle(tri), j in 1:3], color=sol.u[j], colorrange=(0.0, 1.0), colormap=:binary)
                lines!(ax, [p2[get_boundary_nodes(tri)]'; p1[get_boundary_nodes(tri)]'], linewidth=3, color=:black)
            end
            j ≠ 1 && hideydecorations!(ax)
            lad = get_leading_edge(sol.u[j], prob, 1 / 2, centroid)
            if !isempty(lad)
                append!(lad, lad[:, 1])
                if !use_concave_hull
                    if !rotate
                        lines!(ax, lad[1, :], lad[2, :], color=:blue, linewidth=6)
                    else
                        lines!(ax, lad[2, :], lad[1, :], color=:blue, linewidth=6)
                    end
                else
                    ch = concave_hull([[x, y] for (x, y) in eachcol(lad)])
                    verts = ch.vertices
                    if !rotate
                        lines!(ax, verts, color=:blue, linewidth=6)
                    else
                        lad1 = first.(verts)
                        lad2 = last.(verts)
                        lines!(ax, lad2, lad1, color=:blue, linewidth=6)
                    end
                end
            end
            if show_text
                if j == 1
                    if !shift_vals
                        text!(ax, [(25.0, 400.0)]; text=[L"D = %$(rpad(round(D, digits = 5), 2, '0'))"], fontsize=48)
                        text!(ax, [(25.0, 340.0)]; text=[L"λ = %$(rpad(round(λ, digits = 5), 2, '0'))"], fontsize=48)
                        text!(ax, [(25.0, 280.0)]; text=[L"u₀ = %$(rpad(round(u₀, digits = 5), 2, '0'))"], fontsize=48)
                    else
                        text!(ax, [(25.0, 650.0)]; text=[L"D = %$(rpad(round(D, digits = 5), 2, '0'))"], fontsize=48)
                        text!(ax, [(25.0, 650.0 - 60.0)]; text=[L"λ = %$(rpad(round(λ, digits = 5), 2, '0'))"], fontsize=48)
                        text!(ax, [(25.0, 650.0 - 120.0)]; text=[L"u₀ = %$(rpad(round(u₀, digits = 5), 2, '0'))"], fontsize=48)
                    end
                end
            end
        end
    end
    return fig
end

"""
    plot_snapshots(prof, snapshot_save_times;
        product=false,
        fix_D=false,
        fix_λ=false,
        fix_u₀=false,
        fixed_D=nothing,
        fixed_λ=nothing,
        fixed_u₀=nothing,
        use_concave_hull=false,
        rotate=false,
        shift_vals=false)

Plot the snapshots. See also [`get_Dλu₀_for_snapshots`](@ref).

# Arguments
- `prof`: The `ProfileLikelihoodSolution`. See [`find_profiles`](@ref).
- `snapshot_save_times`: The times at which to save the snapshots.

# Keyword Arguments
- `product=false`: Whether to consider `Dλ` or `(D, λ)`.
- `fix_D=false`: Whether to fix `D`.
- `fix_λ=false`: Whether to fix `λ`.
- `fix_u₀=false`: Whether to fix `u₀`.
- `fixed_D=nothing`: The fixed value for `D`, or `nothing` if `fix_D==false`.
- `fixed_λ=nothing`: The fixed value for `λ`, or `nothing` if `fix_λ==false`.
- `fixed_u₀=nothing`: The fixed value for `u₀`, or `nothing` if `fix_u₀==false`.
- `use_concave_hull=false`: Whether to use the concave hull of the leading edge for plotting.
- `rotate=false`: Whether to rotate the mesh by 90 degrees.
- `shift_vals=false`: Whether to shift the values to the top of the figure when showing the values with text.

# Outputs
- `fig`: The figure.
- `_D`: The values of `D` to use for the snapshots, given by `(lower, mle, upper)`. See [`get_Dλu₀_for_snapshots`](@ref).
- `_λ`: The values of `λ` to use for the snapshots, given by `(lower, mle, upper)`. See [`get_Dλu₀_for_snapshots`](@ref).
- `_u₀`: The values of `u₀` to use for the snapshots, given by `(lower, mle, upper)`. See [`get_Dλu₀_for_snapshots`](@ref).
"""
function plot_snapshots(prof, snapshot_save_times;
    product=false,
    fix_D=false,
    fix_λ=false,
    fix_u₀=false,
    fixed_D=nothing,
    fixed_λ=nothing,
    fixed_u₀=nothing,
    use_concave_hull=false,
    rotate=false,
    shift_vals=false)
    mesh = prof.likelihood_problem.data.mesh
    mle_sol = prof.likelihood_solution
    centroid = prof.likelihood_problem.data.centroid
    _D, _λ, _u₀ = get_Dλu₀_for_snapshots(prof, mle_sol;
        product,
        fix_D,
        fix_λ,
        fix_u₀,
        fixed_D,
        fixed_λ,
        fixed_u₀)
    fig = plot_snapshots(_D, _λ, _u₀, mesh, snapshot_save_times, centroid;
        use_concave_hull,
        rotate,
        shift_vals)
    return fig, _D, _λ, _u₀
end

"""
    get_all_leading_edges(_D, _λ, _u₀, mesh, snapshot_save_times, centroid;
        use_concave_hull=false, rotate=false)

Gets the leading edges of the solution at each triplet in `(_D, _λ, _u₀)`.

# Arguments
- `_D`: The `D` values to use. See also [`get_Dλu₀_for_snapshots`](@ref).
- `_λ`: The `λ` values to use. See also [`get_Dλu₀_for_snapshots`](@ref).
- `_u₀`: The `u₀` values to use. See also [`get_Dλu₀_for_snapshots`](@ref).
- `mesh`: The mesh.
- `snapshot_save_times`: The times at which to save the snapshots.
- `centroid`: The centroid of the mesh.

# Keyword Arguments
- `use_concave_hull=false`: Whether to give the concave hull of the leading edge.
- `rotate=false`: Whether to rotate the mesh by 90 degrees.
"""
function get_all_leading_edges(_D, _λ, _u₀, mesh, snapshot_save_times, centroid;
    use_concave_hull=false, rotate=false)
    initial_time = minimum(snapshot_save_times)
    final_time = maximum(snapshot_save_times)
    lads1 = [Vector{Union{Nothing,Vector{Float64}}}(undef, length(snapshot_save_times)) for _ in 1:3]
    lads2 = [Vector{Union{Nothing,Vector{Float64}}}(undef, length(snapshot_save_times)) for _ in 1:3]
    sols = [Any[] for _ in 1:3]
    for i in 1:3
        D = _D[i]
        λ = _λ[i]
        u₀ = _u₀[i]
        BC = get_boundary_conditions(mesh, λ)
        prob, jac_prototype = get_pde(mesh, BC, D, λ, u₀, initial_time, final_time)
        sol = solve(prob,
            TRBDF2(linsolve=KLUFactorization(; reuse_symbolic=true));
            jac_prototype=float.(jac_prototype),
            parallel=true,
            saveat=snapshot_save_times)
        for j in eachindex(snapshot_save_times)
            push!(sols[i], sol.u[j])
            lad = get_leading_edge(sol.u[j], prob, 1 / 2, centroid)
            if isempty(lad)
                lad1, lad2 = nothing, nothing
                lads1[i][j] = lad1
                lads2[i][j] = lad2
            else
                append!(lad, lad[:, 1])
                if !use_concave_hull
                    if !rotate
                        lad1, lad2 = lad[1, :], lad[2, :]
                        lads1[i][j] = lad1
                        lads2[i][j] = lad2
                    else
                        lad1, lad2 = lad[2, :], lad[1, :]
                        lads1[i][j] = lad1
                        lads2[i][j] = lad2
                    end
                else
                    ch = concave_hull([[x, y] for (x, y) in eachcol(lad)])
                    verts = ch.vertices
                    if !rotate
                        lad1, lad2 = first.(verts), last.(verts)
                        lads1[i][j] = lad1
                        lads2[i][j] = lad2
                    else
                        lad1, lad2 = last.(verts), first.(verts)
                        lads1[i][j] = lad1
                        lads2[i][j] = lad2
                    end
                end
            end
        end
    end
    return lads1, lads2, sols
end