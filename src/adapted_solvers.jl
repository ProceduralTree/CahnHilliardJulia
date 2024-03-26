struct adapted_multi_solver <: solver
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

struct adapted_relaxed_multi_solver <: solver
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
