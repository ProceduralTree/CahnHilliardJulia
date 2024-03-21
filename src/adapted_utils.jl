struct adapted_multi_solver
    phase::Matrix{Float64}
    potential::Matrix{Float64}
    xi::Matrix{Float64}
    psi::Matrix{Float64}
    epsilon::Float64
    h::Float64
    dt::Float64
    W_prime::Function
    len::Int
    width::Int

end

function adapted_testgrid(M, len)
    grid = Array{adapted_multi_solver}(undef, len)
    phase = zeros(size(M) .+ 2)
    phase[2:end-1, 2:end-1] = M
    W_prime(x) = -x * (1 - x^2)
    h0 = 3e-3

    for i = 1:len
        grid[i] = adapted_multi_solver(zeros(size(M) .÷ i .+ 2),
            zeros(size(M) .÷ i .+ 2),
            zeros(size(M) .÷ i .+ 2),
            zeros(size(M) .÷ i .+ 2),
            8e-3, h0 * 2^i, 1e-3,
            W_prime,
            size(M, 1) ÷ i, size(M, 2) ÷ i
            )

    end
    copyto!(grid[1].phase, phase)
    return grid
end
