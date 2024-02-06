using ProgressBars
"""
small grid version

Returns
---------------
1 if index i,j is in bounds(without padding) and 0 else
"""
function G(i, j, len, width)
    if 2 <= i <= len + 1 && 2 <= j <= width + 1
        1.0
    else
        0.0
    end
end

function neighbours_in_domain(i, j, len, width)
    (
        G(i + 0.5, j, len, width) +
        G(i - 0.5, j, len, width) +
        G(i, j + 0.5, len, width) +
        G(i, j - 0.5, len, width)
    )

end

"""
discrete laplace operator weighted by boundry to ensure no flux boundry
"""
function discrete_G_weigted_neigbour_sum(i, j, arr, G, len, width)
    (
        G(i + 0.5, j, len, width) * arr[i+1, j] +
        G(i - 0.5, j, len, width) * arr[i-1, j] +
        G(i, j + 0.5, len, width) * arr[i, j+1] +
        G(i, j - 0.5, len, width) * arr[i, j-1]
    )
end

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
