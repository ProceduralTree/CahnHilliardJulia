using JLD2
using LinearAlgebra
using ProgressBars
include("utils.jl")
include("solvers.jl")
include("multisolver.jl")
include("testgrids.jl")
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
