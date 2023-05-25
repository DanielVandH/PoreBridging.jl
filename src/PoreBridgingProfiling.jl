module PoreBridgingProfiling

using DelimitedFiles
using OrdinaryDiffEq
using LinearSolve
using FiniteVolumeMethod
using DelaunayTriangulation
using ProfileLikelihood
using CairoMakie
using ElasticArrays
using Random
using PolygonOps
using LinearAlgebra
using StatsBase
using OptimizationOptimJL
using Optimization
using StaticArrays
using OptimizationNLopt
using LatinHypercubeSampling
using DataFrames
using ConcaveHull
using CSV
const FVM = FiniteVolumeMethod
const DT = DelaunayTriangulation

include("structs.jl")
include("profile_analysis_functions.jl")
include("full_likelihood_analysis_functions.jl")
include("data_setup.jl")
include("geometry_setup.jl")
include("mle.jl")
include("objective_function.jl")
include("optimisation_setup.jl")
include("pde_setup.jl")
include("prediction_intervals.jl")
include("profile.jl")
include("snapshots.jl")
include("summary_statistics.jl")

DelaunayTriangulation.number_type(::Tuple{A,B,C}) where {T,A<:AbstractArray{T},B<:AbstractArray{T},C<:AbstractArray{T}} = T

export get_data_structs,
    OptimisationOptions,
    ProfileAnalysisLikelihoodOptions,
    profile_analysis,
    wave_profile_analysis,
    full_likelihood_analysis,
    wave_full_likelihood_analysis,
    ProfileAnalysisResults,
    ReverseProfileAnalysisResults,
    get_pore_data,
    plot_square_likelihood_analysis,
    plot_wave_likelihood_analysis,
    get_all_leading_edges,
    plot_tissue_growth_predictions,
    get_square_pore_geometry,
    get_wave_pore_geometry,
    get_cross_pore_geometry,
    get_pore_statistics,
    construct_prediction_function,
    pore_prediction_function

end # module