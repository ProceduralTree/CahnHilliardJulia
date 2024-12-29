using DataFrames
using JLD2
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

tests = [testgrid(relaxed_multi_solver , M , 2 , dt = t ) for t in 1e-2./(1:64)]

function iter(g::Vector{T} , n) where T<: solver
    out = []
    for j in 1:n
    set_xi_and_psi!(g[1])
    for i = 1:64
        elyps_solver!(g[1] , 1000)
        alt_v_cycle!(g, 1)
    end
    end
    push!(out, (phase=copy(g[1].phase), iteration=n))
    return out
end


tasks = []
for i in eachindex(tests)
    t = Threads.@spawn iter(tests[i], i)
    push!(tasks , (iteration = 1 , task = t))
    end
result = DataFrame()
for task in tasks
    append!(result , fetch(task.task) )
    end
jldsave("experiments/relaxed-time.jld2"; result)
