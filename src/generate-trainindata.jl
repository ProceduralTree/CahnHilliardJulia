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
SIZES = [2^i for i = 5:9]
SIZE = SIZES[2]
M = fill(-1, SIZE, SIZE)
incirc(M) = filter(x -> norm(x.I .- (size(M, 1) / 2, size(M, 2) / 2)) < min(size(M)...) / 5, CartesianIndices(M))
insquare(M) = filter(x -> norm(x.I .- (size(M, 1) / 2, size(M, 2) / 2), Inf) < min(size(M)...) / 4, CartesianIndices(M))
side(M) = filter(x -> x < CartesianIndex(size(M, 1), size(M, 2) ÷ 2), CartesianIndices(M))
halfcirc(M) = filter(x -> norm(x.I .- (1, size(M, 2) / 2), 2) < 12, CartesianIndices(M))
M[insquare(M)] .= 1
display(heatmap(M, aspect_ratio=:equal, xlims=(1, SIZE)))
t = 50
M = testdata(64, 8, 8, 2)

trainings_input = [
    [testdata(i, 8, 8, j) for i = 2 .^ (5:9), j in [2, Inf]];
    [do M=fill(-1, s, s)
        M[incirc(M)] .= 1
        end  for s = 2 .^ (5:9)
    ]
    [do M=fill(-1, s, s)
        M[insquare(M)] .= 1
        end  for s = 2 .^ (5:9)
    ]
    [do M=fill(-1, s, s)
        M[side(M)] .= 1
        end  for s = 2 .^ (5:9)
    ]
    [do M=fill(-1, s, s)
        M[halfcirc(M)] .= 1
        end  for s = 2 .^ (5:9)
    ]
]


jldsave("data/test-phasefield.jld2"; M=g[1].phase[2:end-1, 2:end-1])

function gen_cycle_data(M::Array)
    out = []
    iter = 2 .^ (0:7)
    g = testgrid(multi_solver, M, 2)
    set_xi_and_psi!(g[1])
    for i = 1:128
        v_cycle!(g, 1)
        if i in iter
            push!(out, (phase=g[1].phase, iteration=i))
        end
    end
    return out
end
