include("solvers.jl")
include("adapted_solvers.jl")
include("utils.jl")
include("multisolver.jl")
include("multi_relaxed.jl")
include("elypssolver.jl")
include("testgrids.jl")
using Plots
using LinearAlgebra
using ProgressBars
using JLD2
using Distributed
M = jldopen("../data/test-phasefield.jld2")["M"]

original_grid = testgrid(multi_solver, M, 2)
alphas = 1e4:1e3:1e6


set_xi_and_psi!(original_grid[1])
for j in 1:10
    v_cycle!(original_grid, 1)
end
tasks = Set{Task}()
for α in alphas
    t = Threads.@spawn alpha_error(α , original_grid[1].phase)
    push!(tasks , t)
end
results  = [fetch(t) for t in tasks]

function alpha_error(α::Number , solution::Array )
    test_solver  = testgrid(relaxed_multi_solver, M, 2, alpha=α)
    set_xi_and_psi!(test_solver[1])
    elyps_solver!(test_solver[1], 1000)
    for j in 1:10
        v_cycle!(test_solver , 1)
    end
return norm(test_solver[1].phase - solution)
end
