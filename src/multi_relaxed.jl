function L(solver::relaxed_multi_solver,i,j , phi , mu)
    xi = solver.phase[i, j] / solver.dt -
         (discrete_G_weigted_neigbour_sum(i, j, solver.potential, G, solver.len, solver.width)
          -
          neighbours_in_domain(i, j, G, solver.len, solver.width) * mu )/solver.h^2
    psi = solver.epsilon^2 * solver.alpha*(solver.c[i,j] - phi) - solver.potential[i,j] - 2 * solver.phase[i,j]
    return [xi, psi]
end

function dL(solver::relaxed_multi_solver , i , j)
    return [ (1/solver.dt) (1/solver.h^2*neighbours_in_domain(i,j,G,solver.len , solver.width));
             (-1*solver.epsilon^2 * solver.alpha  - 2) 1]
    end

function SMOOTH!(
    solver::T,
    iterations,
    adaptive
) where T <: Union{relaxed_multi_solver , adapted_relaxed_multi_solver}
    for k = 1:iterations
        old_phase = copy(solver.phase)
        for I in CartesianIndices(solver.phase)[2:end-1, 2:end-1]
            i, j = I.I
            bordernumber = neighbours_in_domain(i, j, G, solver.len, solver.width)


            solver.phase[I] = (solver.epsilon^2 * solver.alpha * solver.c[I] - solver.h^2 / bordernumber * ( -solver.xi[I]  - discrete_G_weigted_neigbour_sum(i,j,solver.potential , G , solver.len , solver.width) / solver.h^2 ) - solver.psi[I]) / (solver.epsilon^2 * solver.alpha  + 2 + solver.h^2 / (bordernumber*solver.dt))

            #since the solver still needs the potetential we calculate it as well
            solver.potential[I] = (solver.phase[I]/solver.dt - solver.xi[I] - discrete_G_weigted_neigbour_sum(i,j, solver.potential , G , solver.len , solver.width)/solver.h^2) * (-solver.h^2/bordernumber)
        end

        if adaptive && LinearAlgebra.norm(old_phase - solver.phase) < 1e-10
            #println("SMOOTH terminated at $(k) succesfully")
            break
        end
    end
end
