"""
    get_square_pore_geometry(pore_width; max_area_factor=1e-4, min_angle=33.3, max_points=10_000)

Returns the mesh for the square pore geometry.
"""
function get_square_pore_geometry(pore_width; max_area_factor=1e-4, min_angle=33.3, max_points=10_000)
    L = pore_width
    a, b, c, d = 0.0, L, 0.0, L
    points = [(a, c), (b, c), (b, d), (a, d), (a, c)]
    boundary_nodes, points = convert_boundary_points_to_indices(points; existing_points=ElasticMatrix{Float64}(undef, 2, 0))
    tri = triangulate(points; boundary_nodes)
    A = get_total_area(tri)
    refine!(tri, max_area=max_area_factor * A, min_angle=min_angle, max_points=max_points)
    delete_ghost_triangles!(tri)
    mesh = FVMGeometry(tri)
    return mesh
end

"""
    get_wave_pore_geometry(; max_area_factor=1e-4, min_angle=33.3, max_points=10_000)

Returns the mesh for the wave pore geometry.
"""
function get_wave_pore_geometry(; max_area_factor=1e-4, min_angle=33.3, max_points=10_000)
    wave_pore_bnd = readdlm("data/boundaries/PORE_SWV500DAY14IMAGE15PORE1.dat")
    wave_pore_bnd = reverse(wave_pore_bnd; dims=1)
    wave_pore_bnd = unique(wave_pore_bnd; dims=1)
    wave_pore_bnd[:, 1] .-= minimum(wave_pore_bnd[:, 1])
    wave_pore_bnd[:, 2] .-= minimum(wave_pore_bnd[:, 2])
    xmin, xmax = extrema(wave_pore_bnd[:, 1])
    ymin, ymax = extrema(wave_pore_bnd[:, 2])
    c = ((xmin + xmax) / 2, (ymin + ymax) / 2)
    all_angles = [atan(y - c[2], x - c[1]) for (x, y) in eachrow(wave_pore_bnd)]
    sortidx = sortperm(all_angles)
    pt1idx = findfirst(sortidx .== 1)
    circshift!(sortidx, -pt1idx + 1)
    sortidx[[77, 76, 75, 74, 73, 72, 71, 70]] .= sortidx[[75, 70, 72, 74, 73, 76, 77, 71]]
    deleteat!(sortidx, [74, 76, 77])
    wave_pore_bnd = wave_pore_bnd[sortidx, :]
    x = wave_pore_bnd[:, 1]
    y = wave_pore_bnd[:, 2]
    push!(x, x[begin])
    push!(y, y[begin])
    boundary_nodes, points = convert_boundary_points_to_indices(x, y; existing_points=ElasticMatrix{Float64}(undef, 2, 0))
    tri = triangulate(points; boundary_nodes)
    A = get_total_area(tri)
    refine!(tri, max_area=max_area_factor * A, min_angle=min_angle, max_points=max_points)
    delete_ghost_triangles!(tri)
    mesh = FVMGeometry(tri)
    return mesh
end

"""
    get_cross_pore_geometry(L, cut; max_area_factor=1e-4, min_angle=33.3, max_points=10_000)

Returns the mesh for the cross pore geometry.
"""
function get_cross_pore_geometry(L, cut; max_area_factor=1e-4, min_angle=33.3, max_points=10_000)
    A = (cut, 0.0)
    B = (L - cut, 0.0)
    C = (0.0, cut)
    D = (0.0, L - cut)
    E = (cut, cut)
    F = (cut, L - cut)
    G = (L - cut, cut)
    H = (L - cut, L - cut)
    I = (L, cut)
    J = (L, L - cut)
    K = (cut, L)
    M = (L - cut, L)
    pts = [A, B, G, I, J, H, M, K, F, D, C, E, A]
    x = first.(pts)
    y = last.(pts)
    boundary_nodes, points = convert_boundary_points_to_indices(x, y; existing_points=ElasticMatrix{Float64}(undef, 2, 0))
    tri = triangulate(points; boundary_nodes)
    A = get_total_area(tri)
    refine!(tri, max_area=max_area_factor * A, min_angle=min_angle, max_points=max_points)
    delete_ghost_triangles!(tri)
    mesh = FVMGeometry(tri)
    return mesh
end