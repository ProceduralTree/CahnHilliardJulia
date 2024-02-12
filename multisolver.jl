include("/home/jon/Projects/julia_tst/utils.jl")

function L(solver::multi_solver,i,j)
    xi = solver.phase[i, j] / solver.dt -
         (discrete_G_weigted_neigbour_sum(i, j, solver.potential, G, solver.len, solver.width)
          -
          neighbours_in_domain(i, j, G, solver.len, solver.width) * solver.potential[i, j])/solver.h^2
    psi = solver.epsilon^2/solver.h^2 *
          (discrete_G_weigted_neigbour_sum(i, j, solver.phase, G, solver.len, solver.width)
           -
           neighbours_in_domain(i, j, G, solver.len, solver.width) * solver.phase[i, j]) - 2 * solver.phase[i, j] + solver.potential[i, j]
    return [xi, psi]
end

function dL(solver::multi_solver , i , j)
    return [ 1/solver.dt 1/solver.h^2*neighbours_in_domain(i,j,G,solver.len , solver.width);
             -1*solver.epsilon^2/solver.h^2 * (neighbours_in_domain(i,j,G,solver.len , solver.width) - 2) 1]
    end

function SMOOTH!(
    solver::multi_solver,
    iterations,
    adaptive
)
    for k = 1:iterations
        old_phase = copy(solver.phase)
        for I in CartesianIndices(solver.phase)[2:end-1, 2:end-1]
            i, j = I.I
            bordernumber = neighbours_in_domain(i, j, G, solver.len, solver.width)
            coefmatrix =
                [
                    (1/solver.dt) (bordernumber/solver.h^2);
                    (-1*(2+(solver.epsilon^2/solver.h^2)*bordernumber)) 1
                ]


            b =
                [
                    (
                        solver.xi[i, j]
                        +
                        discrete_G_weigted_neigbour_sum(
                            i, j, solver.potential, G, solver.len, solver.width
                        )
                        /
                        solver.h^2
                    ),
                    (
                        solver.psi[i, j]
                        -
                        (solver.epsilon^2 / solver.h^2)
                        *
                        discrete_G_weigted_neigbour_sum(
                            i, j, solver.phase, G, solver.len, solver.width
                        )
                    )
                ]

            res = coefmatrix \ b
            solver.phase[i, j] = res[1]
            solver.potential[i, j] = res[2]

        end

        if adaptive && LinearAlgebra.norm(old_phase - solver.phase) < 1e-8
            print("SMOOTH terminated at $(k) succesfully")
            break
        end
    end
end
