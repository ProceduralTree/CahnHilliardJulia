function SMOOTH(
    xi,
    psi,
    phase,
    mu,
    epsilon,
    h,
    dt,
    len,
    width,
    iterations,
    adaptive
)
    for k =  ProgressBar(1:iterations)
        old_phase = copy(phase)
        for i = 2:(len + 1)
            for j = 2:(width + 1)
                bordernumber = neighbours_in_domain(i, j, len, width)
                coefmatrix =
                    [
                        (1 / dt)  (bordernumber / h^2) ;
                        (-1 * (2 + (epsilon^2 / h^2) * bordernumber))  1
                    ]


                b =
                    [
                        (
                            xi[i, j]
                            + discrete_G_weigted_neigbour_sum(
                                i, j, mu, G, len, width
                            )
                            / h^2
                        ),
                        (
                            psi[i, j]
                            - (epsilon^2 / h^2)
                            * discrete_G_weigted_neigbour_sum(
                                i, j, phase, G, len, width
                            )
                        )
                    ]

                res = coefmatrix \ b
                phase[i, j] = res[1]
                mu[i, j] = res[2]

            end
        end

        if adaptive && LinearAlgebra.norm(old_phase - phase) < 1e-8
            print("SMOOTH terminated at $(k) succesfully")
            break
        end
    end
     (phase, mu)
    end
