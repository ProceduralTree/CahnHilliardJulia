include("solvers.jl")
include("adapted_solvers.jl")
include("utils.jl")
include("multisolver.jl")
include("multi_relaxed.jl")
include("elypssolver.jl")
include("testgrids.jl")
using Plots
using JLD2
using LinearAlgebra
using ProgressBars
using Random
using DataFrames

incirc(M) = filter(x -> norm(x.I .- (size(M, 1) / 2, size(M, 2) / 2)) < min(size(M)...) / 5, CartesianIndices(M))
insquare(M) = filter(x -> norm(x.I .- (size(M, 1) / 2, size(M, 2) / 2), Inf) < min(size(M)...) / 4, CartesianIndices(M))
side(M) = filter(x -> x.I[2] < size(M, 2) ÷ 2, CartesianIndices(M))
halfcirc(M) = filter(x -> norm(x.I .- (1, size(M, 2) / 2), 2) < min(size(M)...) / 5, CartesianIndices(M))

function get_special_input(fn, size)
    M = fill(-1, size , size )
    M[fn(M)] .= 1
    return M
end

function gen_cycle_data(M::Array)
    out = []
    iter = 2 .^ (0:6)
    g = testgrid(multi_solver, M, 2)
    set_xi_and_psi!(g[1])
    for i = 1:128
        v_cycle!(g, 1)
        if i in iter
            push!(out, (phase=copy(g[1].phase), iteration=i))
        end
    end
    out = [merge(o,(;expected=g[1].phase)) for o in out]
    return out
end


t1= [testdata(i, i ÷ 4, i /4 , j) for i = 2 .^ (5:9), j in [2, Inf]]
t2 = [get_special_input(fn,s) for s = 2 .^ (5:9), fn in [halfcirc , incirc, side , insquare]]
trainings_input = reshape([t1;; t2], *(size([t1;; t2])...))
tasks = []
for i in eachindex(trainings_input)
    t = Threads.@spawn gen_cycle_data(trainings_input[i])
    push!(tasks , (iteration = 1 , task = t))
    end
result = DataFrame()
for task in tasks
    append!(result , fetch(task) )
    end
jldsave("data/trainings_data.jld2"; result)
