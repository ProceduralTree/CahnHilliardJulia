include(pwd() * "/src/" * "solvers.jl")
include(pwd() * "/src/" * "utils.jl")
function testgrid(::Type{multi_solver},M, len)
    grid = Array{multi_solver}(undef, len)
    phase = zeros(size(M) .+ 2)
    phase[2:end-1, 2:end-1] = M
    W_prime(x) = -x * (1 - x^2)
    h0 = 3e-3


    for i = 1:len
        grid[i] = multi_solver(zeros(size(M) .÷ i .+ 2),
            zeros(size(M) .÷ i .+ 2),
            zeros(size(M) .÷ i .+ 2),
            zeros(size(M) .÷ i .+ 2),
            8e-3, h0 * 2^i, 1e-3,
            W_prime,
            size(M, 1) ÷ i, size(M, 2) ÷ i)

    end
    copyto!(grid[1].phase, phase)
    return grid

end

function testgrid(::Type{relaxed_multi_solver},M, len)
    grid = Array{relaxed_multi_solver}(undef, len)
    phase = zeros(size(M) .+ 2)
    phase[2:end-1, 2:end-1] = M
    W_prime(x) = -x * (1 - x^2)
    h0 = 3e-3

    for i = 1:len
        grid[i] = relaxed_multi_solver(zeros(size(M) .÷ i .+ 2),
            zeros(size(M) .÷ i .+ 2),
            zeros(size(M) .÷ i .+ 2),
            zeros(size(M) .÷ i .+ 2),
            zeros(size(M) .÷ i .+ 2),
            8e-3, h0 * 2^i, 1e-3,
            W_prime,
            size(M, 1) ÷ i, size(M, 2) ÷ i,
            1000001)

    end
    copyto!(grid[1].phase, phase)
    return grid
end
