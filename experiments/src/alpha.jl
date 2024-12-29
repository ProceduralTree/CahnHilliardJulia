include(pwd() * "/src/solvers.jl")
include(pwd() * "/src/adapted_solvers.jl")
include(pwd() * "/src/utils.jl")
include(pwd() * "/src/multisolver.jl")
include(pwd() * "/src/multi_relaxed.jl")
include(pwd() * "/src/testgrids.jl")
include(pwd() * "/src/elypssolver.jl")
using Plots
using LaTeXStrings
using LinearAlgebra
using Printf
using ProgressBars
default(fontfamily="computer modern" , titlefontsize=32 , guidefontsize=22 , tickfontsize = 22 , legendfontsize=22)
pgfplotsx()
layout2x2 = grid(2,2)
layout3x1 = @layout [ b  c ; a]
size3x1 = (1600,1600)
SIZE = 64
M = testdata(SIZE, SIZE ÷ 5, SIZE /5 , 2)

using JLD2
using Distributed
using ProgressBars
using DataFrames

original_grid = testgrid(multi_solver, M, 2)
alphas = 0:1e4:2e6

function alpha_error(alpha::Number , solution::Array )
    test_solver  = testgrid(relaxed_multi_solver, M, 2, alpha=alpha)
    set_xi_and_psi!(test_solver[1])
    for j in 1:64
        elyps_solver!(test_solver[1], 1000)
        alt_v_cycle!(test_solver , 1)
    end
return [(;alpha=alpha , error=norm(test_solver[1].phase - solution))]
end
set_xi_and_psi!(original_grid[1])
for j in 1:64
    alt_v_cycle!(original_grid, 1)
end
print("finished original v_cycle")
tasks = []
for alpha in alphas
    t = Threads.@spawn alpha_error(alpha , original_grid[1].phase)
    push!(tasks , (alpha=alpha , task = t))
end
result = DataFrame()
for task in ProgressBar(tasks)
    append!(result , fetch(task.task) )
    end
jldsave("experiments/alpha.jld2"; result)
