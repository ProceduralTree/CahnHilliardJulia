struct multi_solver
    phase::Matrix{Float64}
    potential::Matrix{Float64}
    xi::Matrix{Float64}
    psi::Matrix{Float64}
    epsilon::Float64
    h::Float64
    dt::Float64
    W_prime::Function
    len::UInt32
    width::UInt32

end

function testdata(gridsize , blobs , radius ,norm)
rngpoints = rand(1:gridsize, 2, blobs)
M = zeros(gridsize,gridsize) .- 1
for p in axes(rngpoints , 2)
    point = rngpoints[:, p]
    for I in CartesianIndices(M)
        if (LinearAlgebra.norm(point .- I.I  , norm) < radius)
            M[I] = 1
        end
    end
end
M

end

function set_xi_and_psi!(solver::multi_solver)
    xi_init(x) = x / solver.dt
    psi_init(x) = solver.W_prime(x) - 2 * x
    solver.xi[2:end-1, 2:end-1] = xi_init.(solver.phase[2:end-1,2:end-1])
    solver.psi[2:end-1, 2:end-1] = psi_init.(solver.phase[2:end-1,2:end-1])
    return nothing
end
