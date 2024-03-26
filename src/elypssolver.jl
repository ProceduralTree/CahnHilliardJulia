using ProgressBars

"""
    elyps_solver(c,
    phase,
    len,
    width,
    alpha,
    h,
    n
)

TBW
"""

function elyps_solver!(solver::T, n) where T  <: Union{relaxed_multi_solver , adapted_relaxed_multi_solver}
    for k in 1:n
        for i = 2:(solver.len+1)
            for j = 2:(solver.width+1)
                bordernumber = neighbours_in_domain(i, j,G, solver.len, solver.width)
                solver.c[i, j] =
                    (
                        solver.alpha * solver.phase[i, j] +
                        discrete_G_weigted_neigbour_sum(i, j, solver.c, G, solver.len, solver.width) / solver.h^2
                    ) / (bordernumber / solver.h^2 + solver.alpha)

            end
        end
    end
end
