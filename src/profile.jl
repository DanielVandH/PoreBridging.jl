"""
    find_profiles(likprob, mle_sol, resolution=50, min_steps=10, maxiters=100)

Computes the profile likelihoods.

# Arguments
- `likprob`: The `LikelihoodProblem`.
- `mle_sol`: The maximum likelihood solution. See [`find_mle`](@ref).
- `resolution`: The number of points to use in each direction when computing the profile likelihoods.
- `min_steps`: The minimum number of steps in each direction to use when computing the profile likelihoods.
- `maxiters`: The maximum number of iterations in the optimiser to use when computing the profile likelihoods.

# Outputs
- `prof`: The `ProfileLikelihoodSolution`.
"""
function find_profiles(likprob, mle_sol, resolution=50, min_steps=10, maxiters=100)
    @time prof = profile(likprob,
        mle_sol;
        ftol_abs=1e-8,
        ftol_rel=1e-8,
        xtol_abs=1e-8,
        xtol_rel=1e-8,
        maxiters,
        maxtime=600,
        parallel=true,
        resolution=resolution,
        next_initial_estimate_method=:interp,
        min_steps=min_steps)
    return prof
end
