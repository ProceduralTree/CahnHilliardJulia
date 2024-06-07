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
default(fontfamily="computer modern" , titlefontsize=32 , guidefontsize=32 , tickfontsize = 22 )
pgfplotsx()
layout2x2 = grid(2,2)
layout3x1 = @layout [ b  c ; a]
size3x1 = (1600,1600)
SIZE = 64
M = testdata(SIZE, SIZE ÷ 5, SIZE /5 , 2)

using JLD2
using DataFrames
results = jldopen("experiments/iteration.jld2")["result"]
anim = @animate for res in eachrow(results)
    heatmap(res.solver.phase , title="phase field" , legend=:none , aspectratio=:equal , showaxis=false , grid=false , size=(400 ,400))
end
gif(anim , "images/iteration.gif" , fps = 10)
