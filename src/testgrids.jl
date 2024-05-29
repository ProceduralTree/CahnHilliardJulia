function W_prime(x)
    return -x * (1 - x^2)
end
function testgrid(::Type{multi_solver},M, len; dt = 1e-3 ,  epsilon=8e-3 , h0=3e-3)
    grid = Array{multi_solver}(undef, len)
    phase = zeros(size(M) .+ 2)
    phase[2:end-1, 2:end-1] = M


    for i = 1:len
        dims = size(M) .÷ 2^(i-1) .+ 2
        grid[i] = multi_solver(zeros(dims),
            zeros(dims),
            zeros(dims),
            zeros(dims),
            epsilon, h0 * 2^i, dt,
            W_prime,
            (dims .- 2)...)

    end
    copyto!(grid[1].phase, phase)
    return grid

end

function testgrid(::Type{relaxed_multi_solver},M, len ; alpha=1e6 , dt=1e-3, epsilon=8e-3 , h0=3e-3)
    grid = Array{relaxed_multi_solver}(undef, len)
    phase = zeros(size(M) .+ 2)
    phase[2:end-1, 2:end-1] = M

    for i = 1:len
        dims = size(M) .÷ 2^(i-1) .+ 2
        grid[i] = relaxed_multi_solver(zeros(dims),
            zeros(dims),
            zeros(dims),
            zeros(dims),
            zeros(dims),
            epsilon, h0 * 2^i, dt,
            W_prime,
            (dims .- 2)... ,
            alpha)

    end
    copyto!(grid[1].phase, phase)
    return grid
end
