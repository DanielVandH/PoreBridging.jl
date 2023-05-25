################################################################################################
##
## SETUP 
##
################################################################################################
using PoreBridgingProfiling
using Setfield
using JLD2
using CairoMakie
using ProfileLikelihood

jld2_path_profiles = normpath(@__DIR__, "..", "paper", "jld2_files", "profile_results")
jld2_path_full_likelihood = normpath(@__DIR__, "..", "paper", "jld2_files", "full_likelihood_results")
jld2_path_figures = normpath(@__DIR__, "..", "paper", "figures", "reparametrisation")

square_perimeter_area_full_likelihood_results_1 = load(joinpath(jld2_path_full_likelihood, "square_perimeter_area_full_likelihood_results_1.jld2"))["square_perimeter_area_full_likelihood_results"];
square_area_perimeter_profile_results = load(joinpath(jld2_path_profile, "square_area_perimeter_results.jld2"))["square_area_perimeter_results"];

threshold = ProfileLikelihood.get_chisq_threshold(0.95, 2)

################################################################################################
##
## GET GRIDS
##
################################################################################################
lb = [1e-6, 1e-6]
ub_product = [500.0, 5.0]
ub_separate = [1000.0, 5.0]
regular_grid_product = RegularGrid(lb, ub_product, 40)
regular_grid_separate = RegularGrid(lb, ub_separate, 40)
grid_1_prod = get_range(regular_grid_product, 1)
grid_2_prod = get_range(regular_grid_product, 2)
grid_1_sep = get_range(regular_grid_separate, 1)
grid_2_sep = get_range(regular_grid_separate, 2)

################################################################################################
##
## PLOTTING
##
################################################################################################
product_vals = square_perimeter_area_full_likelihood_results_1.loglik_vals_product
separate_vals = square_perimeter_area_full_likelihood_results_1.loglik_vals_separate

fig = Figure(fontsize=33, resolution=(1550, 555))
ax = Axis(fig[1, 1],
    xlabel=L"$D$ [$\mu$m$^2$/day]",
    ylabel=L"$\lambda$ [1/day]",
    title=L"(a): $(D, \lambda)$ parametrisation",
    titlealign=:left,
    width=650,
    height=400,
    xticks=(0:100:500, [L"%$s" for s in 0:100:500]),
    yticks=(0:1:5, [L"%$s" for s in 0:1:5]))
contourf!(ax, grid_1_sep, grid_2_sep, separate_vals, colormap=:viridis, levels=-10:0.25:0, extendlow = :auto)
contour!(ax, grid_1_sep, grid_2_sep, separate_vals, levels=[threshold], linewidth=3, color=:red)
xlims!(ax, 0, 500)
ylims!(ax, 0, 5)

ax = Axis(fig[1, 2],
    xlabel=L"$D\lambda$ [$\mu$m$^2$/day$^2]",
    ylabel=L"$\lambda$ [1/day]",
    title=L"(b): $(D\lambda, \lambda)$ parametrisation",
    titlealign=:left,
    width=650,
    height=400,
    xticks=(0:100:500, [L"%$s" for s in 0:100:500]),
    yticks=(0:1:5, [L"%$s" for s in 0:1:5])
)
contourf!(ax, grid_1_prod, grid_2_prod, product_vals, colormap=:viridis, levels=-10:0.25:0, extendlow = :auto)
contour!(ax, grid_1_prod, grid_2_prod, product_vals, levels=[threshold], linewidth=3, color=:red)
xlims!(ax, 0, 500)
ylims!(ax, 0, 5)

save(joinpath(jld2_path_figures, "reparametrisation_contours.pdf"), fig)