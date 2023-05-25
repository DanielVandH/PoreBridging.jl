"""
    get_parameter_bounds(; geo=:square)

Returns the bounds on the parameters `(D, λ, u₀)`. If `geo==:square`,
the returned values are:

- `(Dlow, Dupp, λlow, λupp, u₀low, u₀upp) = (1e-6, 1000.0, 1e-6, 5.0, 1e-6, 0.5)`,

otherwise, if `geo==:wave`, then 

- `(Dlow, Dupp, λlow, λupp, u₀low, u₀upp) = (1e-6, 2000.0, 1e-6, 10.0, 1e-6, 0.5)`.
"""
function get_parameter_bounds(; geo=:square)
    if geo == :square
        Dlow = 1e-6
        Dupp = 2000.0
        λlow = 1e-6
        λupp = 5.0
        u₀low = 1e-6
        u₀upp = 0.5
    else
        Dlow = 1e-6
        Dupp = 2000.0
        λlow = 1e-6
        λupp = 10.0
        u₀low = 1e-6
        u₀upp = 0.5
    end
    return Dlow, Dupp, λlow, λupp, u₀low, u₀upp
end

"""
    construct_likelihood_problem(mesh, initial_time, final_time, saveat_times,
        area_data, perimeter_data, radius_data, sigma,
        pore_area, pore_perimeter, pore_smallest_rad, centroid;
        include_area=true, include_perimeter=true, include_radius=true, product=false,
        fix_D, fix_λ, fix_u₀, fixed_D=nothing, fixed_λ=nothing, fixed_u₀=nothing,
        geo=:square)

Returns the `LikelihoodProblem` using the given options.

# Arguments
- `mesh`: The mesh.
- `initial_time`: The initial time of the data.
- `final_time`: The final time of the data.
- `saveat_times`: The times to save the numerical solutions at for the objective function.
- `area_data`: The void coverages.
- `perimeter_data`: The normalised void perimeters.
- `radius_data`: The normalised void radii.
- `sigma`: The noise estimates.
- `pore_area`: The pore area.
- `pore_perimeter`: The pore perimeter.
- `pore_smallest_rad`: The smallest distance from the `centroid` to the pore boundary.
- `centroid`: The centroid of the pore.

# Keyword Arguments
- `include_area`: Whether to include the void coverages in the log-likelihood function.
- `include_perimeter`: Whether to include the normalised void perimeters in the log-likelihood function.
- `include_radius`: Whether to include the normalised void radii in the log-likelihood function.
- `product`: Whether to consider `Dλ` or `(D, λ)`.
- `fix_D`: Whether to fix `D`.
- `fix_λ`: Whether to fix `λ`.
- `fix_u₀`: Whether to fix `u₀`.
- `fixed_D=nothing`: The fixed value for `D`, or `nothing` if `fix_D==false`.
- `fixed_λ=nothing`: The fixed value for `λ`, or `nothing` if `fix_λ==false`.
- `fixed_u₀=nothing`: The fixed value for `u₀`, or `nothing` if `fix_u₀==false`.
- `geo=:square`: The geometry to use. Either `:square` or `:wave`.
"""
function construct_likelihood_problem(mesh, initial_time, final_time, saveat_times,
    area_data, perimeter_data, radius_data, sigma,
    pore_area, pore_perimeter, pore_smallest_rad, centroid;
    include_area=true, include_perimeter=true, include_radius=true, product=false,
    fix_D, fix_λ, fix_u₀, fixed_D=nothing, fixed_λ=nothing, fixed_u₀=nothing,
    geo=:square)
    fix_D && isnothing(fixed_D) && throw("No value for D provided, but D is fixed.")
    fix_λ && isnothing(fixed_λ) && throw("No value for λ provided, but λ is fixed.")
    fix_u₀ && isnothing(fixed_u₀) && throw("No value for u₀ provided, but u₀ is fixed.")
    ## Get the parameters
    p = (
        mesh=mesh,
        initial_time=initial_time,
        final_time=final_time,
        saveat_times=saveat_times,
        area_data=area_data,
        perimeter_data=perimeter_data,
        radius_data=radius_data,
        sigma=sigma,
        pore_area=pore_area,
        pore_perimeter=pore_perimeter,
        pore_smallest_rad=pore_smallest_rad,
        centroid=centroid,
        include_area=include_area,
        include_perimeter=include_perimeter,
        include_radius=include_radius
    )
    ## Get the bounds 
    Dlow, Dupp, λlow, λupp, u₀low, u₀upp = get_parameter_bounds(; geo)
    lb = [Dlow, λlow, u₀low]
    ub = [Dupp, λupp, u₀upp]
    ## Get the initial estimates 
    θ₀ = (lb .+ ub) ./ 2
    ## Now define the problem based on what we are fixing, and whether we are taking a product
    if !product
        if !fix_D && !fix_λ && !fix_u₀
            θ_conv = θ -> θ
            _lb = lb
            _ub = ub
            _θ₀ = θ₀
            syms = [:D, :λ, :u₀]
        elseif !fix_D && !fix_λ && fix_u₀
            θ_conv = θ -> [θ[1], θ[2], fixed_u₀]
            _lb = lb[[1, 2]]
            _ub = ub[[1, 2]]
            _θ₀ = θ₀[[1, 2]]
            syms = [:D, :λ]
        elseif !fix_D && fix_λ && !fix_u₀
            θ_conv = θ -> [θ[1], fixed_λ, θ[2]]
            _lb = lb[[1, 3]]
            _ub = ub[[1, 3]]
            _θ₀ = θ₀[[1, 3]]
            syms = [:D, :u₀]
        elseif fix_D && !fix_λ && !fix_u₀
            θ_conv = θ -> [fixed_D, θ[1], θ[2]]
            _lb = lb[[2, 3]]
            _ub = ub[[2, 3]]
            _θ₀ = θ₀[[2, 3]]
            syms = [:λ, u₀]
        elseif !fix_D && fix_λ && fix_u₀
            θ_conv = θ -> [θ[1], fixed_λ, fixed_u₀]
            _lb = lb[[1]]
            _ub = ub[[1]]
            _θ₀ = θ₀[[1]]
            syms = [:D]
        elseif fix_D && !fix_λ && fix_u₀
            θ_conv = θ -> [fixed_D, θ[1], fixed_u₀]
            _lb = lb[[2]]
            _ub = ub[[2]]
            _θ₀ = θ₀[[2]]
            syms = [:λ]
        elseif fix_D && fix_λ && !fix_u₀
            θ_conv = θ -> [fixed_D, fixed_λ, θ[1]]
            _lb = lb[[3]]
            _ub = ub[[3]]
            _θ₀ = θ₀[[3]]
            syms = [:u₀]
        elseif fix_D && fix_λ && fix_u₀
            throw("Cannot fix all three parameters.")
        end
    elseif product
        if !fix_D && !fix_λ && !fix_u₀
            θ_conv = θ -> [θ[1] / θ[2], θ[2], θ[3]]
            _lb = [lb[1] * lb[2], lb[2], lb[3]]
            _ub = [ub[1] * ub[2], ub[2], ub[3]]
            _θ₀ = [θ₀[1] * θ₀[2], θ₀[2], θ₀[3]]
            syms = [:Dλ, :λ, :u₀]
        elseif !fix_D && !fix_λ && fix_u₀
            θ_conv = θ -> [θ[1] / θ[2], θ[2], fixed_u₀]
            _lb = [lb[1] * lb[2], lb[2]]
            _ub = [ub[1] * ub[2], ub[2]]
            _θ₀ = [θ₀[1] * θ₀[2], θ₀[2]]
            syms = [:Dλ, :λ]
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
    loglik_fnc = (θ, p) -> pore_likelihood(θ_conv(θ), p)
    likprob = LikelihoodProblem(
        loglik_fnc,
        _θ₀;
        data=p,
        syms=syms,
        f_kwargs=(adtype=Optimization.AutoFiniteDiff(),),
        prob_kwargs=(lb=_lb, ub=_ub)
    )
    return likprob
end