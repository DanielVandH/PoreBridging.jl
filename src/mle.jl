"""
    perform_grid_search(likprob, n=500, gens=1000)

Performs a grid search on the likelihood problem `likprob` using 
Latin hypercube sampling, evaluating the likelihood at `n` parameters set 
from a Latin hypercube with `1000` generations (see `LHCoptim`). Returns 
`(gs_ir, loglik_vals)`, where `gs_ir` are the grid search results, 
and `loglik_vals` are the associated values of the log-likelihood.
"""
function perform_grid_search(likprob, n=500, gens=1000)
    d = length(likprob.syms)
    plan, _ = LHCoptim(n, d, gens; threading=false)
    lb = get_lower_bounds(likprob)
    ub = get_upper_bounds(likprob)
    bnds = [(ℓ, u) for (ℓ, u) in zip(lb, ub)]
    parameter_vals = Matrix(scaleLHC(plan, bnds)')
    irregular_grid = IrregularGrid(lb, ub, parameter_vals)
    @time gs_ir, loglik_vals = grid_search(
        likprob,
        irregular_grid;
        save_vals=Val(true),
        parallel=Val(true)
    )
    return gs_ir, loglik_vals
end

"""
    perform_grid_search(likprob, regular_grid::RegularGrid)

Performs a grid search on the likelihood problem `likprob` using
the regular grid `regular_grid`. Returns `(gs_ir, loglik_vals)`, where
`gs_ir` are the grid search results, and `loglik_vals` are the associated
values of the log-likelihood.
"""
function perform_grid_search(likprob, regular_grid::RegularGrid)
    return grid_search(likprob, regular_grid; save_vals=Val(true),parallel=Val(true))
end

"""
    find_mle(likprob, gs_ir, alg=NLopt.LN_BOBYQA)

Computes the MLEs for the given likelihood problem `likprob`.

# Arguments 
- `likprob`: The `LikelihoodProblem`.
- `gs_ir`: The grid search results. See [`perform_grid_search`](@ref).
- `alg`: The algorithm to use for the optimisation.

# Outputs 
- `refined_likprob`: The refined likelihood problem, using the results from the grid search.
- `mle_sol`: The maximum likelihood solution.
""" 
function find_mle(likprob, gs_ir, alg=NLopt.LN_BOBYQA)
    refined_likprob = update_initial_estimate(likprob, gs_ir)
    @time mle_sol = mle(
        refined_likprob,
        alg;
        ftol_abs=1e-4,
        ftol_rel=1e-4,
        xtol_abs=1e-4,
        xtol_rel=1e-4
    )
    return refined_likprob, mle_sol
end