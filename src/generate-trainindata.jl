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
SIZES = [2^i for i=5:9 ]
SIZE = SIZES[2]
M = fill(-1,SIZE,SIZE)
incirc(M) = filter(x-> norm(x.I .- (size(M,1)/2 , size(M,2)/2)) < 12 , CartesianIndices(M))
insquare(M) = filter(x-> norm(x.I .- (size(M,1)/2 , size(M,2)/2), Inf) < 12 , CartesianIndices(M))
side(M) = filter(x-> x < CartesianIndex(size(M,1),size(M,2) ÷ 2) , CartesianIndices(M))
halfcirc(M) = filter(x-> norm(x.I .- (1 , size(M,2)/2), 2) < 12 , CartesianIndices(M))
M[halfcirc(M)] .= 1
display(heatmap(M , aspect_ratio=:equal , xlims = (1,SIZE)))
t = 50
M = testdata(64, 8, 8, 2)
g = testgrid(M, 2)
pbar = ProgressBar(total = t * 10)
for k = 1:t
    set_xi_and_psi!(g[1])
    for l = 1:100
        v_cycle!(g,1)
        update(pbar)
    end
end
jldsave("data/test-phasefield.jld2";M = g[1].phase[2:end-1,2:end-1])
