################################################################################################
##
## SETUP PACKAGES AND DATA
##
################################################################################################
using PoreBridgingProfiling
using DelaunayTriangulation
using CairoMakie
using DelimitedFiles
using JLD2
using ProfileLikelihood

jld2_path_full_likelihood = normpath(@__DIR__, "..", "paper", "jld2_files", "full_likelihood_results")
jld2_path_profiles = normpath(@__DIR__, "..", "paper", "jld2_files", "profile_results")
jld2_path_figures = normpath(@__DIR__, "..", "paper", "figures", "hypothetical")
jld2_path_hypothetical = normpath(@__DIR__, "..", "paper", "jld2_files", "hypothetical_study")

reference_likelihood_results = load(joinpath(jld2_path_full_likelihood, "square_area_full_likelihood_results_1.jld2"))["square_area_full_likelihood_results"];
square_area_perimeter_profile_results = load(joinpath(jld2_path_profiles, "square_area_perimeter_results.jld2"))["square_area_perimeter_results"];
reference_profile_results = square_area_perimeter_profile_results[1]

############################################
## SETUP THE GEOMETRY
############################################
L = 525.0
cut = L / 3.5
mesh = get_cross_pore_geometry(L, cut)

############################################
## GET THE SUMMARY STATISTICS
############################################
pore_area, pore_perimeter, pore_smallest_rad, centroid = get_pore_statistics(mesh)
snapshot_save_times = [5, 7, 14, 25, 28]

############################################
## PLOT THE SNAPSHOTS
############################################
_D = reference_profile_results.C_D
_λ = reference_profile_results.C_λ
_u₀ = reference_profile_results.C_u₀
lad1, lad2, sols = get_all_leading_edges(_D, _λ, _u₀, mesh, snapshot_save_times, centroid)
tri = mesh.mesh_information.triangulation
pts = tri.points
bn = tri.boundary_nodes
pts_bn = pts[:, bn]
b1 = pts_bn[1, :]
b2 = pts_bn[2, :]
fig = plot_tissue_growth_predictions(lad1, lad2, sols, b1, b2, pts, tri, snapshot_save_times)
save(joinpath(jld2_path_figures, "cross_tissue_growth_snapshots.pdf"), fig)

############################################
## COMPUTE THE SUMMARY STATISTICS
############################################
## Need to start by finding the feasible points from the original likelihood 
Dλ_grid, λ_grid = reference_likelihood_results.param_ranges
loglik_vals = reference_likelihood_results.loglik_vals_product
threshold = ProfileLikelihood.get_chisq_threshold(0.95, 2)
feasible_idx = findall(>(threshold), loglik_vals)
feasible_Dλ = [Dλ_grid[I[1]] for I in feasible_idx]
feasible_λ = [λ_grid[I[2]] for I in feasible_idx]
feasible_pairs = [(Dλ, λ) for (Dλ, λ) in zip(feasible_Dλ, feasible_λ)]

## Now we setup the parameters for propagating the uncertainty 
nt = 100
full_q_vals = zeros(3nt + 1, length(feasible_pairs))
initial_time = 5.0
final_time = 14.0
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

## Now we can compute the prediction function everywhere in the credible region 
if isfile(joinpath(jld2_path_hypothetical, "hypothetical_q_vals.jld2"))
    full_q_vals = load(joinpath(jld2_path_hypothetical, "hypothetical_q_vals.jld2"))["full_q_vals"]
else
    Base.Threads.@threads for (i, θ) in collect(enumerate(feasible_pairs))
        @show i, θ
        full_q_vals[:, i] .= pore_prediction_function([θ..., _u₀[1]], p)
    end
    jldsave(joinpath(jld2_path_hypothetical, "hypothetical_q_vals.jld2"); full_q_vals)
end

## Extract the uncertainty bands 
q_mle = pore_prediction_function([get_mle(reference_likelihood_results.mle_sol_product)..., _u₀[1]], p)
q_lwr = minimum(full_q_vals; dims=2) |> vec
q_upr = maximum(full_q_vals; dims=2) |> vec
intervals = [(l, u) for (l, u) in zip(q_lwr, q_upr)]
coverage_intervals = intervals[1:nt]
perimeter_intervals = intervals[(nt+1):2nt]
bridging_sample = full_q_vals[end, :]

## Now let's also plot the summary statistics
fig = Figure(fontsize=60)
ax1 = Axis(fig[1, 1],
    xlabel=L"$t$ [day]",
    ylabel=L"\mu_{\mathrm{c}}(t)",
    title=L"(a):$ $ Void coverage",
    titlealign=:left,
    width=600,
    height=400,
    xticks=(0:10:40, [L"%$s" for s in 0:10:40]),
    yticks=(0:0.25:1, [L"%$s" for s in 0:0.25:1]))
coverage_ℓ = first.(coverage_intervals)
coverage_u = last.(coverage_intervals)
coverage_mle = q_mle[1:nt]
lines!(ax1, saveat_times_fine, coverage_ℓ, linewidth=4, color=:black)
lines!(ax1, saveat_times_fine, coverage_u, linewidth=4, color=:black)
lines!(ax1, saveat_times_fine, coverage_mle, linewidth=4, color=:black, linestyle=:dash)
xlims!(ax1, 0, 42)
ylims!(ax1, 0, 1)

ax2 = Axis(fig[1, 2],
    xlabel=L"$t$ [day]",
    ylabel=L"\mu_{\mathrm{p}}(t)",
    title=L"(b):$ $ Void perimeter",
    titlealign=:left,
    width=600,
    height=400,
    xticks=(0:10:40, [L"%$s" for s in 0:10:40]),
    yticks=(0:0.25:1, [L"%$s" for s in 0:0.25:1]))
perimeter_ℓ = first.(perimeter_intervals)
perimeter_u = last.(perimeter_intervals)
perimeter_mle = q_mle[(nt+1):2nt]
lines!(ax2, saveat_times_fine, perimeter_ℓ, linewidth=4, color=:black)
lines!(ax2, saveat_times_fine, perimeter_u, linewidth=4, color=:black)
lines!(ax2, saveat_times_fine, perimeter_mle, linewidth=4, color=:black, linestyle=:dash)
xlims!(ax2, 0, 42)
ylims!(ax2, 0, 1)

ax3 = Axis(fig[1, 3],
    xlabel=L"$t_b$ [day]",
    ylabel=L"p(t_b)",
    title=L"(c): $p(t_b)$",
    titlealign=:left,
    width=600,
    height=400,
    xticks=(11:3:29, [L"%$s" for s in 11:3:29]),
    yticks=(0:0.1:0.2, [L"%$s" for s in 0:0.1:0.2]))
density!(ax3, bridging_sample, strokecolor=:black, color=(:black, 0.0), strokewidth=4, strokearound=false)
vlines!(ax3, [q_mle[end]], color=:black, linewidth=4, linestyle=:dash)
resize_to_layout!(fig)
save(joinpath(jld2_path_figures, "cross_tissue_growth_summary_statistics.pdf"), fig)