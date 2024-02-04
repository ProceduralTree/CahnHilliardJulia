include("smooth.jl")
include("elypssolver.jl")
using RestrictProlong


"""
    restrict!(smallgrid_solver::multi_solver , largegrid_solver::multi_solver)::multi_solver

------------
Requires
----------
smallgrid solver and largegid solvers to be multiple of 2 from each other bar padding eg. (66x66)->(34x34)


TBW
"""
function restrict!(smallgrid_solver::multi_solver , largegrid_solver::multi_solver)::multi_solver
   largegrid_solver.phase[2:end-1,2:end-1] = restrict(smallgrid_solver.phase[2:end-1,2:end-1])
   largegrid_solver.potential[2:end-1,2:end-1] = restrict(smallgrid_solver.potential[2:end-1,2:end-1])
   largegrid_solver.xi[2:end-1,2:end-1] = restrict(smallgrid_solver.xi[2:end-1,2:end-1])
   largegrid_solver.psi[2:end-1,2:end-1] = restrict(smallgrid_solver.psi[2:end-1,2:end-1])
    largegrid_solver.h = 2 * smallgrid_solver.h
    largegrid_solver.dt =  smallgrid_solver.dt
    largegrid_solver.epsilon = smallgrid_solver.epsilon
    largegrid_solver.W_prime = smallgrid_solver.W_prime
    return nothing
    end
function prolong!(smallgrid::multi_solver , largegrid::multi_solver)

end

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


function L(solver::multi_solver , i , j )
xi = solver.phase[i,j] / solver.dt - (discrete_G_weigted_neigbour_sum(i,j, solver.potential[i,j] , G , solver.len , solver.width)  - neighbours_in_domain(i,j,solver.len , solver.width) * solver.potential[i,j]) / solver.h^2
psi = solver.epsilon^2 *
    (discrete_G_weigted_neigbour_sum(i,j, solver.phase[i,j] , G , solver.len , solver.width) / h^2
     - neighbours_in_domain(i,j,solver.len , solver.width) * solver.phase[i,j]) - 2 * solver.phase[i,j] + solver.potential[i,j]
    return (xi , psi)
end





function v_cycle(grid::Array{multi_solver}, level)

    SMOOTH(400)
    solver = grid[level]

    # extract (d,r) as array operations

    dr = zeros((solver.len+ 2, solver.width+ 2, 2))

    # TODO check array indicies
    dr[2:end-1, 2:end-1, :] = [
            [
               [solver.xi[i,j] , solver.psi[i,j]] .- L(solver ,i,j)
                for j in 2:(solver.width+ 1)
            ]
            for i in 2:(solver.len+ 1)
        ]
    )
    d = dr[:, :, 0]
    r = dr[:, :, 1]

    # print(f"Max derivation d: {np.linalg.norm(d)}")
    # print(f"Max derivation r: {np.linalg.norm(r)}")
    restrict!(grid[level], grid[level + 1])
    solver = grid[level + 1]

    u_large = zeros((solver.len+ 2, solver.width+ 2))
    v_large = zeros((solver.len+ 2, solver.width+ 2))

    # solve for phi^ mu^ with L
    for i in 2:(self.len+1)
        for j in 2:(self.width+1)

            # print(f"Max derivation u: {np.linalg.norm(u_large)}")
            # print(f"Max derivation v: {np.linalg.norm(v_large)}")
            end
        end


# smooth again:
SMOOTH(800)
end
