###############################################################################
#                  Common Utility Functions For Multi Solvers                 #
###############################################################################

"""
restricts an array on the small grid to an array in the large grid asserts size arr=2^n + 2 and returns ret=2^(n-1) + 2

Returns
---------------------------
large grid array + padding
"""
function restrict(arr, G)
    shape = (size(arr) .- 2) .÷ 2
    ret = zeros(shape .+ 2)
    for I in CartesianIndices(ret)[2:end-1, 2:end-1]
        i, j = I.I
        g = [
            G(2 * i - 1, 2 * j - 1, (size(arr) .- 2)...),
            G(2 * i - 1, 2 * j, (size(arr) .- 2)...),
            G(2 * i, 2 * j - 1, (size(arr) .- 2)...),
            G(2 * i, 2 * j, (size(arr) .- 2)...)
        ]
        if sum(g) == 0
            ret[I] = 0
        else
            ret[I] = (
                1 / sum(g)
                *
                dot(g,
                    [
                        arr[2*i-1, 2*j-1],
                        arr[2*i-1, 2*j],
                        arr[2*i, 2*j-1],
                        arr[2*i, 2*j]
                    ]
                )
            )
        end
    end
    return ret
end

"""
    prolong(arr , G)

interpolates int a smaller grid by a factor of 2

"""
function prolong(arr, G)
    inner_shape = (size(arr) .- 2) .* 2
    ret = zeros(inner_shape .+ 2)
    ONE = oneunit(CartesianIndices(arr)[1])
    for I in CartesianIndices(arr)[2:end-1, 2:end-1]
        Ind = 2 * (I - ONE) + ONE
        for J in (Ind-ONE):Ind
            ret[J] = G(J.I..., inner_shape...) * arr[I]
        end
    end
    return ret
end

struct multi_solver
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
struct relaxed_multi_solver
    phase::Matrix{Float64}
    potential::Matrix{Float64}
    xi::Matrix{Float64}
    psi::Matrix{Float64}
    c::Matrix{Float64}
    epsilon::Float64
    h::Float64
    dt::Float64
    W_prime::Function
    len::Int
    width::Int
    alpha::Float64

end

function testgrid(M, len)
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

function relaxed_testgrid(M, len)
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

"""
    restrict!(smallgrid_solver::multi_solver , largegrid_solver::multi_solver)::multi_solver

------------
Requires
----------
smallgrid solver and largegid solvers to be multiple of 2 from each other bar padding eg. (66x66)->(34x34)

------------
Returns
------------
    nothing. mutatest largegid in place to represent the smallgrid

"""
function restrict_solver!(smallgrid_solver::T, largegrid_solver::T) where {T<:Union{multi_solver,relaxed_multi_solver}}
    copy!(largegrid_solver.phase, restrict(smallgrid_solver.phase, G))
    copy!(largegrid_solver.potential, restrict(smallgrid_solver.potential, G))
    return nothing
end

"""
    neighbours_in_domain(i, j, G, len, width)

TBW
counts neigbours in domain
"""
function neighbours_in_domain(i, j, G, len, width)
    (
        G(i + 0.5, j, len, width)
        + G(i - 0.5, j, len, width)
        + G(i, j + 0.5, len, width)
        + G(i, j - 0.5, len, width)
    )

end

"""
    discrete_G_weigted_neigbour_sum(i, j, arr, G, len, width)

TBW
--------
sums all neighbours depending on wheter tey are in the domain determined by G
"""
function discrete_G_weigted_neigbour_sum(i, j, arr, G, len, width)
    (
        G(i + 0.5, j, len, width) * arr[i+1, j]
        + G(i - 0.5, j, len, width) * arr[i-1, j]
        + G(i, j + 0.5, len, width) * arr[i, j+1]
        + G(i, j - 0.5, len, width) * arr[i, j-1]
    )
end

function testdata(gridsize , blobs , radius ,norm)
rngpoints = rand(1:gridsize, 2, blobs)
M = zeros(gridsize,gridsize) .- 1
for p in axes(rngpoints , 2)
    point = rngpoints[:, p]
    for I in CartesianIndices(M)
        if (LinearAlgebra.norm(point .- I.I  , norm) < radius)
            M[I] = 1
        end
    end
end
M

end

"""
Boundry indicator function

Returns
---------------
1 if index i,j is in bounds(without padding) and 0 else
"""
function G(i, j, len, width)
    if 2 <= i <= len + 1 && 2 <= j <= width + 1
        return 1.0
    else
        return 0.0
    end
end

function set_xi_and_psi!(solver::T) where T <: Union{multi_solver , relaxed_multi_solver}
    xi_init(x) = x / solver.dt
    psi_init(x) = solver.W_prime(x) - 2 * x
    solver.xi[2:end-1, 2:end-1] = xi_init.(solver.phase[2:end-1,2:end-1])
    solver.psi[2:end-1, 2:end-1] = psi_init.(solver.phase[2:end-1,2:end-1])
    return nothing
end

function bulk_energy(solver::T) where T <: Union{multi_solver , relaxed_multi_solver}
    energy = 0
    dx = CartesianIndex(1,0)
    dy = CartesianIndex(0,1)
    W(x) = 1/4 * (1-x^2)^2
    for I in CartesianIndices(solver.phase)[2:end-1,2:end-1]
        i,j = I.I
        energy += solver.epsilon^2 / 2 * G(i+ 0.5,j ,solver.len, solver.width) * (solver.phase[I+dx] - solver.phase[I])^2 + G(i,j+0.5,solver.len ,solver.width) * (solver.phase[I+dy] - solver.phase[I])^2 + W(solver.phase[I])
        end
   return energy
end
