using DataFrames
using JLD2
using ProgressMeter
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


M = testdata(2^10 , 2^5 , 2^7 , 2 )
grids = testgrid(relaxed_multi_solver  , M , 7 , h0=3e-3 * 64 /1024)
# inits
for i=2:size(grids,1)
    restrict_solver!(grids[i-1] , grids[i])
end
tests = [[grids[i-1] , grids[i]] for i=2:size(grids,1)]

n = 4
m = 1024

function iter(g::Vector{T} , n , prg::Progress) where T<: solver
    out = []
    for j in 1:n
    set_xi_and_psi!(g[1])
    elyps_solver!(g[1] , 1000)
    for i = 1:m
        alt_v_cycle!(g, 1)
        next!(prg)
    end
    push!(out, (phase=copy(g[1].phase), iteration=j))
    end
    return out
end


prg=Progress(size(tests ,1)*n*m , showspeed=true , )
tasks = []
for i in eachindex(tests)
    t = Threads.@spawn iter(tests[i], n , prg)
    push!(tasks , (iteration = 1 , task = t))
    end
result = DataFrame()
for task in tasks
    append!(result , fetch(task.task) )
    end
jldsave("experiments/relaxed_space_refinement.jld2"; result)
