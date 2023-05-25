"""
    pore_likelihood(area_data, perimeter_data, radius_data,
        area_cache, perimeter_cache, radius_cache::AbstractArray{T},
        include_area, include_perimeter, include_radius,
        sigma) where {T}

Returns the log-likelihood function for the given data. 

# Arguments 
- `area_data`: The void coverages.
- `perimeter_data`: The normalised void perimeters.
- `radius_data`: The normalised void radii.
- `area_cache`: The void coverages for the numerical solutions.
- `perimeter_cache`: The normalised void perimeters for the numerical solutions.
- `radius_cache`: The normalised void radii for the numerical solutions.
- `include_area`: Whether to include the void coverages in the log-likelihood function.
- `include_perimeter`: Whether to include the normalised void perimeters in the log-likelihood function.
- `include_radius`: Whether to include the normalised void radii in the log-likelihood function.
- `sigma`: The noise estimates.

# Outputs 
- `ℓ`: The log-likelihood function for the given data.
"""
function pore_likelihood(area_data, perimeter_data, radius_data,
    area_cache, perimeter_cache, radius_cache::AbstractArray{T},
    include_area, include_perimeter, include_radius,
    sigma) where {T}
    ℓ = zero(T)
    for j in eachindex(area_data)
        area_sigma = get_σ(sigma, j, :area)
        perimeter_sigma = get_σ(sigma, j, :perimeter)
        radius_sigma = get_σ(sigma, j, :radius)
        μ_area = area_cache[j]
        μ_perimeter = perimeter_cache[j]
        μ_radius = radius_cache[j]
        areas = @views area_data[j]
        perimeters = @views perimeter_data[j]
        radii = @views radius_data[j]
        for i in eachindex(area_data[j])
            if include_area
                ℓ = ℓ - log(area_sigma * sqrt(2π)) - 0.5 * (areas[i] - μ_area)^2 / area_sigma^2
            end
            if include_perimeter
                ℓ = ℓ - log(perimeter_sigma * sqrt(2π)) - 0.5 * (perimeters[i] - μ_perimeter)^2 / perimeter_sigma^2
            end
            if include_radius
                ℓ = ℓ - log(radius_sigma * sqrt(2π)) - 0.5 * (radii[i] - μ_radius)^2 / radius_sigma^2
            end
        end
    end
    return ℓ
end

"""
    pore_likelihood(θ::AbstractVector{T}, param) where {T}

Returns the log-likelihood function for the given parameters `θ = (D, λ, u₀)` and parameters 
`param` from the `LikelihoodProblem`.
"""
function pore_likelihood(θ::AbstractVector{T}, param) where {T}
    D, λ, u₀ = θ
    (; mesh, initial_time, final_time, saveat_times,
        area_data, perimeter_data, radius_data, sigma,
        pore_area, pore_perimeter, pore_smallest_rad, centroid,
        include_area, include_perimeter, include_radius) = param
    BC = get_boundary_conditions(mesh, λ)
    prob, jac_prototype = get_pde(mesh, BC, D, λ, u₀, initial_time, final_time)
    sol = solve(prob,
        TRBDF2(linsolve=KLUFactorization(; reuse_symbolic=false));
        jac_prototype=float.(jac_prototype),
        saveat=saveat_times,
        parallel=true)
    if !SciMLBase.successful_retcode(sol)
        return typemin(T)
    end
    area_cache = zeros(length(saveat_times))
    perimeter_cache = zeros(length(saveat_times))
    radius_cache = zeros(length(saveat_times))
    threshold_parameter = 1 / 2
    compute_summary_statistics!(area_cache, perimeter_cache, radius_cache,
        sol, threshold_parameter, prob, centroid,
        pore_area, pore_perimeter, pore_smallest_rad)
    ℓ = pore_likelihood(area_data, perimeter_data, radius_data,
        area_cache, perimeter_cache, radius_cache,
        include_area, include_perimeter, include_radius, sigma)
    @show ℓ, θ
    return ℓ
end

"""
    get_σ(σ, time_idx, type)

Returns the estimate for the standard deviation from `σ`
for the given summary statistic `type` at the time 
index `time_idx`.
"""
get_σ(σ::Number, time_idx, type) = σ
get_σ(σ, time_idx, type) = σ[type][time_idx]