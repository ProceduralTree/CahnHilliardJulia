
using LinearAlgebra
using Plots
include("utils.jl")
include("smooth.jl")
include("multisolver.jl")
include("elypssolver.jl")


M = testdata(512, 20, 40, 2)

display(heatmap(t))

I = CartesianIndex(2, 2)
phase = zeros(size(M) .+ 2);
phase[2:end-1, 2:end-1] = M;
mu = copy(phase);
W_prime(x) = -x * (1 - x^2)
solver = multi_solver(
    phase,
    zeros(size(phase)),
    zeros(size(phase)),
    zeros(size(phase)),
    8e-3, 1e-3, 1e-3,
    W_prime,
    size(M, 1), size(M, 2))
set_xi_and_psi!(solver)
SMOOTH!(solver, 400, true)

d = zeros(size(solver.phase))
r = zeros(size(solver.phase))

for I in CartesianIndices(solver.phase)[2:end-1, 2:end-1]
    d[I], r[I] = [solver.xi[I], solver.psi[I]] .- L(solver, I.I...)
end

display(heatmap(d))

function v_cycle(grid::Array{multi_solver}, level)

    SMOOTH!(solver, 400, true)
    solver = grid[level]

    # extract (d,r) as array operations

    d = zeros(size(solver.phase))
    r = zeros(size(solver.phase))

    for I in CartesianIndices(solver.phase)[2:end-1, 2:end-1]
        d[I], r[I] = [solver.xi[I], solver.psi[I]] .- L(solver, I.I...)
    end

    # print(f"Max derivation d: {np.linalg.norm(d)}")
    # print(f"Max derivation r: {np.linalg.norm(r)}")
    restrict!(grid[level], grid[level+1])
    solver = grid[level+1]

    u_large = zeros((solver.len + 2, solver.width + 2))
    v_large = zeros((solver.len + 2, solver.width + 2))
    #TODO short newton iteration for
    for i = 1:3
        for I in CartesianIndices(solver.phase)[2:end-1, 2:end-1]

        end
    end

    u_large = zeros((solver.len + 2, solver.width + 2))
    v_large = zeros((solver.len + 2, solver.width + 2))

    # solve for phi^ mu^ with L
    for i in 2:(solver.len+1)
        for j in 2:(solver.width+1)

            # print(f"Max derivation u: {np.linalg.norm(u_large)}")
            # print(f"Max derivation v: {np.linalg.norm(v_large)}")
        end
    end
    # smooth again:
    SMOOTH!(solver, 800, true)
end
