include("elypssolver.jl")
function SMOOTH!(
    solver::multi_solver,
    iterations,
    adaptive
)
    for k = 1:iterations
        old_phase = copy(solver.phase)
        for i = 2:(solver.len + 1)
            for j = 2:(solver.width + 1)
                bordernumber = neighbours_in_domain(i, j, solver.len, solver.width)
                coefmatrix =
                    [
                        (1 / solver.dt)  (bordernumber / solver.h^2) ;
                        (-1 * (2 + (solver.epsilon^2 / solver.h^2) * bordernumber))  1
                    ]


                b =
                    [
                        (
                            solver.xi[i, j]
                            + discrete_G_weigted_neigbour_sum(
                                i, j, solver.potential, G, solver.len, solver.width
                            )
                            / solver.h^2
                        ),
                        (
                            solver.psi[i, j]
                            - (solver.epsilon^2 / solver.h^2)
                            * discrete_G_weigted_neigbour_sum(
                                i, j, solver.phase, G, solver.len, solver.width
                            )
                        )
                    ]

                res = coefmatrix \ b
                solver.phase[i, j] = res[1]
                solver.potential[i, j] = res[2]

            end
        end

        if adaptive && LinearAlgebra.norm(old_phase - solver.phase) < 1e-3
            print("SMOOTH terminated at $(k) succesfully")
            break
        end
    end
    end
