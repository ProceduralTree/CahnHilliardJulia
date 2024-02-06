function L(solver::multi_solver , i , j )
xi = solver.phase[i,j] / solver.dt -
    (discrete_G_weigted_neigbour_sum(i,j, solver.potential , G , solver.len , solver.width)  - neighbours_in_domain(i,j,solver.len , solver.width) * solver.potential[i,j]) / solver.h^2
psi = solver.epsilon^2 *
    (discrete_G_weigted_neigbour_sum(i,j, solver.phase , G , solver.len , solver.width) / solver.h^2
     - neighbours_in_domain(i,j,solver.len , solver.width) * solver.phase[i,j]) - 2 * solver.phase[i,j] + solver.potential[i,j]
    return (xi , psi)
end

function dL(solver::multi_solver , i , j)
    return [ 1/solver.dt neighbours_in_domain(i,j,G);
             (neighbours_in_domain(i,j,G) - 2) 1]
    end
