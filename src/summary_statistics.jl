"""
    get_smallest_distance(xy::AbstractArray{T}, centroid::AbstractArray{T})

Given a matrix of coordinates `xy` defining some boundary, with each column representing 
a coordinate on the boundary, and the coordinates are in counter-clockwise order, 
returns the smallest distance from `centroid` to the boundary curve.
"""
function get_smallest_distance(xy::AbstractArray{T}, centroid::AbstractArray{T}) where {T}
    d = typemax(T)
    N = size(xy, 2)
    for i in 2:N
        @views new_d = smallest_distance(xy[:, i-1], xy[:, i], centroid)
        d = new_d < d ? new_d : d
    end
    @views new_d = smallest_distance(xy[:, end], xy[:, 1], centroid)
    d = new_d < d ? new_d : d
    return d
end

"""
    smallest_distance(p₁, p₂, q)

Given a line segment `p₁p₂` and a point `q`, returns the smallest distance from `q` to the line segment.

https://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment
"""
@inline function smallest_distance(p₁::AbstractArray{T}, p₂::AbstractArray{T}, q) where {T} # smallest distance from q to the line segment p₁₂
    qp₁ = zero(p₁)
    p₁p₂ = zero(p₂)
    return smallest_distance!(qp₁, p₁p₂, p₁, p₂, q)
end
@inline function smallest_distance!(qp₁::AbstractArray{T}, p₁p₂, p₁, p₂, q) where {T}
    @. qp₁ = q - p₁
    @. p₁p₂ = p₂ - p₁
    t = dot(qp₁, p₁p₂) / dot(p₁p₂, p₁p₂)
    ts = min(max(t, zero(T)), one(T))
    @. p₁p₂ = p₁ + ts * p₁p₂
    @. qp₁ = q - p₁p₂
    return norm(qp₁)
end

"""
    compute_summary_statistics!(areas, perimeters, radii, sol, threshold, prob, centroid, pore_area, pore_perimeter, pore_smallest_rad)

Computes the void coverages, normalised void perimeters, and normalised smallest radii of the void.

# Arguments 
- `areas`: Cache vector for the void coverages. The new void coverages will be placed into here. 
- `perimeters`: Cache vector for the normalised void perimeters. The new perimeters will be placed into here. 
- `radii`: Cache vector for the smallest radii. The new smallest radii will be placed into here. 
- `saveat_times`: Times to compute the statistics at. 
- `sol`: A solution to the PDE. 
- `threshold`: The threshold parameter `τ` that defines the leading edge. 
- `prob`: The `FVMProblem` defining the PDE. 
- `centroid`: Centroid of the void.
- `pore_area`: Area of the pore. 
- `pore_perimeter`: Perimeter of the pore. 
- `pore_smallest_rad`: Smallest radii of the pore. 

# Outputs 
There are no outputs - all results are placed into `areas`, `perimeters`, and `radii`.
"""
function compute_summary_statistics!(areas, perimeters, radii, sol, threshold, prob, centroid, pore_area, pore_perimeter, pore_smallest_rad)
    for i in eachindex(sol)
        A, P, R = @views compute_leading_edge_statistics(sol.u[i], FiniteVolumeMethod.get_elements(prob), FiniteVolumeMethod.get_points(prob), threshold, centroid)
        areas[i] = one(eltype(sol)) - A / pore_area
        perimeters[i] = P == zero(eltype(sol)) ? one(eltype(sol)) : P / pore_perimeter
        radii[i] = isinf(R) ? one(eltype(sol)) : R / pore_smallest_rad
    end
    return nothing
end

"""
    compute_leading_edge_statistics(u, T, pts, threshold, centroid)

Given a vector of weights `u` for elements `T` with coordinates `pts`, a threshold 
`threshold`, and a reference centroid `centroid`, computes the area of the part of 
the solution above `threshold`, the length of the boundary curve, and the smallest 
distance from `centroid` to the boundary curve.

It is assumed that the leading edge is simply connected.

See also [`get_leading_edge`](@ref).
"""
function compute_leading_edge_statistics(u, T, pts, threshold, centroid)
    A = zero(eltype(u))
    ℓ = zero(eltype(u))
    rad = typemax(eltype(u))
    p = zeros(eltype(u), 2)
    q = zeros(eltype(u), 2)
    r = zeros(eltype(u), 2)
    s = zeros(eltype(u), 2)
    for T in T
        i, j, k = DelaunayTriangulation.indices(T)
        pᵢ = pts[:, i]
        pⱼ = pts[:, j]
        pₖ = pts[:, k]
        uᵢ, uⱼ, uₖ = u[i], u[j], u[k]
        _a, _ℓ, rad = threshold_intersection_statistics!(p, q, r, s, rad, centroid, threshold, uᵢ, uⱼ, uₖ, pᵢ, pⱼ, pₖ)
        A += _a
        ℓ += _ℓ
    end
    return A, ℓ, rad
end

"""
    get_leading_edge(u, prob, threshold, centroid)

Extracts the leading edge, returning it in counter-clockwise order. The leading edge is represented 
as a matrix, with the `x`-coordinates in the first row and the `y`-coordinates in the second row.

See also [`compute_leading_edge_statistics`](@ref).
"""
function get_leading_edge(u, prob, threshold, centroid)
    ## Get the leading edge 
    leading_edge = ElasticMatrix{Float64}(undef, 2, 0)
    sizehint!(leading_edge, (2, num_points(prob)))
    for (v₁, v₂) in get_edges(FiniteVolumeMethod.get_triangulation(prob))
        if !DelaunayTriangulation.is_ghost_edge(v₁, v₂)
            u₁, u₂ = u[v₁], u[v₂]
            if threshold_intersection_exists(threshold, u₁, u₂)
                p₁ = get_points(prob)[:, v₁]
                p₂ = get_points(prob)[:, v₂]
                t = threshold_intersection(threshold, u₁, u₂)
                intersection_point = @. p₁ + t * (p₂ - p₁)
                append!(leading_edge, intersection_point)
            end
        end
    end
    ## Sort the leading edge 
    θ = zeros(size(leading_edge, 2))
    for j in axes(leading_edge, 2)
        x, y = leading_edge[:, j]
        θ[j] = atan(y - gety(centroid), x - getx(centroid))
    end
    sort_idx = sortperm(θ)
    Base.permutecols!!(leading_edge, sort_idx)
    return leading_edge
end

"""
    threshold_intersection_statistics!(p, q, r, s, d, c, τ::T, u₁, u₂, u₃, p₁, p₂, p₃) where {T}

Given a triangle `(p₁, p₂, p₃)`, corresponding weights `(u₁, u₂, u₃)`, and a threshold 
`τ`, computes the area of the part of the lifted triangle (the one in the `z`-plane with `z`-coordinates 
given by the weights) that is above the threshold plane, and also the length of the intersection curve. 
The vectors `p`, `q`, `r`, and `s` are cache vectors that are the same type and size as the points `(p₁, p₂, p₃)`.

We also have an argument `d` that is assumed to be the current smallest distance from a query point c. Where applicable, 
we compute the distance to the edge found from the intersection above to `c` and, if it is smaller than this, we return this 
new `d`, otherwise `d` is returned as-is.
"""
function threshold_intersection_statistics!(p, q, r, s, d, c, τ::T, u₁, u₂, u₃, p₁, p₂, p₃) where {T}
    triangle_pos_relative_to_threshold_plane = threshold_intersection_exists(τ, u₁, u₂, u₃)
    if triangle_pos_relative_to_threshold_plane == -1
        return zero(T), zero(T), d
    elseif triangle_pos_relative_to_threshold_plane == 1
        return PolygonOps.area((p₁, p₂, p₃, p₁)) |> abs, zero(T), d
    else
        if (u₁ ≤ τ) && (u₂ ≤ τ) && (u₃ ≥ τ)
            t₂₃ = threshold_intersection(τ, u₂, u₃)
            t₃₁ = threshold_intersection(τ, u₃, u₁)
            eval_line_segment!(p, t₂₃, p₂, p₃)
            eval_line_segment!(q, t₃₁, p₃, p₁)
            A = PolygonOps.area((q, p, p₃, q))
        elseif (u₁ ≤ τ) && (u₂ ≥ τ) && (u₃ ≤ τ)
            t₁₂ = threshold_intersection(τ, u₁, u₂)
            t₂₃ = threshold_intersection(τ, u₂, u₃)
            eval_line_segment!(p, t₁₂, p₁, p₂)
            eval_line_segment!(q, t₂₃, p₂, p₃)
            A = PolygonOps.area((p, p₂, q, p))
        elseif (u₁ ≤ τ) && (u₂ ≥ τ) && (u₃ ≥ τ)
            t₁₂ = threshold_intersection(τ, u₁, u₂)
            t₃₁ = threshold_intersection(τ, u₃, u₁)
            eval_line_segment!(p, t₁₂, p₁, p₂)
            eval_line_segment!(q, t₃₁, p₃, p₁)
            A = abs(PolygonOps.area((p₁, p₂, p₃, p₁))) - abs(PolygonOps.area((p, q, p₁, p)))
        elseif (u₁ ≥ τ) && (u₂ ≤ τ) && (u₃ ≤ τ)
            t₁₂ = threshold_intersection(τ, u₁, u₂)
            t₃₁ = threshold_intersection(τ, u₃, u₁)
            eval_line_segment!(p, t₁₂, p₁, p₂)
            eval_line_segment!(q, t₃₁, p₃, p₁)
            A = PolygonOps.area((q, p₁, p, q))
        elseif (u₁ ≥ τ) && (u₂ ≤ τ) && (u₃ ≥ τ)
            t₁₂ = threshold_intersection(τ, u₁, u₂)
            t₂₃ = threshold_intersection(τ, u₂, u₃)
            eval_line_segment!(p, t₁₂, p₁, p₂)
            eval_line_segment!(q, t₂₃, p₂, p₃)
            A = abs(PolygonOps.area((p₁, p₂, p₃, p₁))) - abs(PolygonOps.area((q, p, p₂, q)))
        elseif (u₁ ≥ τ) && (u₂ ≥ τ) && (u₃ ≤ τ)
            t₂₃ = threshold_intersection(τ, u₂, u₃)
            t₃₁ = threshold_intersection(τ, u₃, u₁)
            eval_line_segment!(p, t₂₃, p₂, p₃)
            eval_line_segment!(q, t₃₁, p₃, p₁)
            A = abs(PolygonOps.area((p₁, p₂, p₃, p₁))) - abs(PolygonOps.area((p, p₃, q, p)))
        else
            @show u₁, u₂, u₃
        end
    end
    d′ = smallest_distance!(r, s, p, q, c)
    new_d = d′ < d ? d′ : d
    @. p = p - q
    ℓ = norm(p)
    return abs(A), ℓ, new_d
end

"""
    threshold_intersection_exists(τ, u₁, u₂, u₃)

Given a threshold `τ` and a triangle with nodal values `u₁`, `u₂`, and `u₃`, returns 

- `-1`: If all the `uᵢ` are below `τ`.
- `1`: If all the `uᵢ` are above `τ`.
- `0`: If an intersection between one of the triangle's edges with the threshold plane `u = τ` exists.
"""
function threshold_intersection_exists(τ, u₁, u₂, u₃)
    if (u₁ < τ) && (u₂ < τ) && (u₃ < τ)
        return -1 # below
    elseif (u₁ > τ) && (u₂ > τ) && (u₃ > τ)
        return 1 # above 
    else
        return 0
    end
end

"""
    threshold_intersection(τ, uᵢ, uⱼ)

Returns the value `t` such that `uᵢ + t(uⱼ - uᵢ) = τ`.
"""
@inline threshold_intersection(τ, uᵢ, uⱼ) = (τ - uᵢ) / (uⱼ - uᵢ)

"""
    eval_line_segment!(p, t, pᵢ, pⱼ) 

Evaluates the line segment `pᵢ → pⱼ` at `t`, storing the new coordinate in `p`, i.e. 
puts into `p` the value `pᵢ + t(pⱼ - pᵢ)`.
"""
@inline eval_line_segment!(p, t, pᵢ, pⱼ) = @. p = pᵢ + (pⱼ - pᵢ) .* t

"""
    threshold_intersection_exists(τ, uᵢ, uⱼ)

Returns `true` if there exists a `t ∈ [0, 1]` such that `uᵢ + t(uⱼ - uᵢ) = τ`.
"""
@inline threshold_intersection_exists(τ, uᵢ, uⱼ) = (uᵢ < τ && uⱼ > τ) || (uᵢ > τ && uⱼ < τ)
