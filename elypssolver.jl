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
function elyps_solver(c, phase, len, width, alpha, h, n)
    for k in ProgressBar(1:n)
        for i = 2:(len+1)
            for j = 2:(width+1)
                bordernumber = neighbours_in_domain(i, j, len, width)
                c[i, j] =
                    (
                        alpha * phase[i, j] +
                        discrete_G_weigted_neigbour_sum(i, j, c, G, len, width) / h^2
                    ) / (bordernumber / h^2 + alpha)

            end
        end
    end
    c
end
