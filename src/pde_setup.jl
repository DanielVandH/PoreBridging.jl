"""
    get_boundary_conditions(mesh, λ)

Returns the boundary conditions for the given `mesh` and and proliferation rate `λ`.
"""
function get_boundary_conditions(mesh, λ)
    boundary_function = (x, y, t, u, p) -> p.λ * u * (1 - u)
    type = :dudt
    params = [(λ=λ,)]
    BC = BoundaryConditions(mesh, boundary_function, type; params)
    return BC
end

"""
    get_pde(mesh, BC, D, λ, u₀, initial_time, final_time)

Returns `(prob, jac)`, where `prob` is the the PDE problem for the given `mesh`, 
boundary conditions `BC`, diffusivity `D`, proliferation rate `λ`, initial condition `u₀`, 
`initial_time`, and `final_time`. `jac` is the prototype for the Jacobians for the PDE solver.
"""
function get_pde(mesh, BC, D, λ, u₀, initial_time, final_time)
    reaction_function = (x, y, t, u, p) -> p.λ * u * (1.0 - u)
    function flux_function!(q, x, y, t, α, β, γ, p)
        u = α * x + β * y + γ
        q[1] = -p.D * u * α
        q[2] = -p.D * u * β
        return nothing
    end
    flux_parameters = (D=D,)
    reaction_parameters = (λ=λ,)
    initial_condition = zeros(num_points(mesh.mesh_information.triangulation))
    initial_condition[get_boundary_nodes(mesh.mesh_information.triangulation)] .= u₀
    prob = FVMProblem(
        mesh,
        BC;
        flux_function=flux_function!,
        reaction_function=reaction_function,
        flux_parameters=flux_parameters,
        reaction_parameters=reaction_parameters,
        initial_time=initial_time,
        final_time=final_time,
        initial_condition=initial_condition
    )
    jac_prototype = jacobian_sparsity(prob)
    return prob, jac_prototype
end