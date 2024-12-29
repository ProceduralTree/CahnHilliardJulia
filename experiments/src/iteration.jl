using JLD2
using DataFrames
using Random
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

incirc(M) = filter(x -> norm(x.I .- (size(M, 1) / 2, size(M, 2) / 2)) < min(size(M)...) / 3, CartesianIndices(M))
insquare(M) = filter(x -> norm(x.I .- (size(M, 1) / 2, size(M, 2) / 2), Inf) < min(size(M)...) / 4, CartesianIndices(M))
side(M) = filter(x -> x.I[2] < size(M, 2) ÷ 2, CartesianIndices(M))
halfcirc(M) = filter(x -> norm(x.I .- (1, size(M, 2) / 2), 2) < min(size(M)...) / 3, CartesianIndices(M))

function get_special_input(fn, size)
    M = fill(-1, size , size )
    M[fn(M)] .= 1
    return M
end
SIZE  =64
t1= [testdata(SIZE, SIZE ÷ 5, SIZE /5 , j) for j in [1,2, Inf]]
t2 = [get_special_input(fn,SIZE) for  fn in [halfcirc , incirc, side , insquare]]
initial_data = [t1 ; t2]
tests = [testgrid(multi_solver, M , 2) for M in initial_data]

function iter(g::Vector{T} , n) where T<: solver
    out = []
    for j in 1:64
    set_xi_and_psi!(g[1])
    for i = 1:64
        v_cycle!(g, 1)
    end
    push!(out, (solver=deepcopy(g[1]), iteration=j , experiment=n))
    end
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
jldsave("experiments/iteration.jld2"; result)
