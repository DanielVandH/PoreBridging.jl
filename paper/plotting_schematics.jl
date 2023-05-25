################################################################################################
##
## SETUP PACKAGES
##
################################################################################################
using PoreBridgingProfiling
using DelaunayTriangulation
using CairoMakie
using DelimitedFiles
jld2_path_figures = normpath(@__DIR__, "..", "paper", "figures", "schematics")

################################################################################################
##
## PLOT EACH SCHEMATIC
##
################################################################################################
pore_width = 475.0
L = 525.0
cut = L / 3.5
mesh = get_square_pore_geometry(pore_width; max_area_factor=0.001, min_angle=30.0, max_points=3000);
wave_mesh = get_wave_pore_geometry(max_area_factor=0.001, min_angle=30.0, max_points=3000);
cross_mesh = get_cross_pore_geometry(L, cut, max_area_factor=0.001, min_angle=30.0, max_points=3000);

fig = Figure(fontsize=48)
ax = Axis(fig[1, 1],
    width=900, height=900,
    title=L"(a):$ $ Square geometry",
    titlealign=:left,
    titlesize=66)
tri = mesh.mesh_information.triangulation
triplot!(ax, tri; show_convex_hull=false, strokewidth=1.8)
lines!(ax, tri.points[1, get_boundary_nodes(tri)], tri.points[2, get_boundary_nodes(tri)], color=:red, linewidth=5)
xlims!(ax, -50, pore_width + 50)
ylims!(ax, -50, pore_width + 50)
text!(ax, 200, pore_width + 1; text=L"\partial\Omega", fontsize=60)
poly!(ax, [(60.0, 60.0), (100.0, 60.0), (100.0, 100.0), (60.0, 100.0), (60.0, 60.0)], color=:white, strokecolor=:black, strokewidth=0.6, overdraw=true)
text!(ax, 66.0, 61.0; text=L"\Omega", fontsize=60)
hidedecorations!(ax)

ax = Axis(fig[1, 2],
    width=900, height=900,
    title=L"(b):$ $ Wave geometry",
    titlealign=:left,
    titlesize=66)
tri = wave_mesh.mesh_information.triangulation
p1 = deepcopy(tri.points[1, :])
p2 = deepcopy(tri.points[2, :])
tri.points[1, :] .= p2
tri.points[2, :] .= p1
triplot!(ax, tri; show_convex_hull=false, strokewidth=1.8)
text!(ax, 25, pore_width + 10; text=L"\partial\Omega", fontsize=60)
lines!(ax, tri.points[1, get_boundary_nodes(tri)], tri.points[2, get_boundary_nodes(tri)], color=:red, linewidth=5)
poly!(ax, [(335.0, 335.0), (390.0, 335.0), (390.0, 390.0), (335.0, 390.0), (335.0, 335.0)], color=:white, strokecolor=:black, strokewidth=0.6, overdraw=true)
text!(ax, 343.0, 336.0; text=L"\Omega", fontsize=60)
hidedecorations!(ax)

ax = Axis(fig[1, 3],
    width=900, height=900,
    title=L"(c):$ $ Cross geometry",
    titlealign=:left,
    titlesize=66)
tri = cross_mesh.mesh_information.triangulation
triplot!(ax, tri; show_convex_hull=false, strokewidth=1.8)
text!(ax, 25, L - cut + 10; text=L"\partial\Omega", fontsize=60)
lines!(ax, tri.points[1, get_boundary_nodes(tri)], tri.points[2, get_boundary_nodes(tri)], color=:red, linewidth=5)
poly!(ax, [(2cut - 10, 2cut - 10), (2cut + 40, 2cut - 10), (2cut + 40, 2cut + 40), (2cut - 10, 2cut + 40), (2cut - 10, 2cut - 10)], color=:white, strokecolor=:black, strokewidth=0.6, overdraw=true)
text!(ax, 2cut + 2, 2cut - 2; text=L"\Omega", fontsize=60)
hidedecorations!(ax)

resize_to_layout!(fig)

fig

save(joinpath(jld2_path_figures, "schematic.png"), fig)