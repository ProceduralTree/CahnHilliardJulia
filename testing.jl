
using Plots
using LinearAlgebra


include("elypssolver.jl")
include("smooth.jl")

SIZE = 2^9
radius = 50
rngpoints = rand(1:SIZE, 2, 10)
M = zeros(SIZE,SIZE) .- 1


for p in axes(rngpoints , 2)
    point = rngpoints[:,p]
    for I in eachindex(IndexCartesian(), M)
            if (LinearAlgebra.norm(point .- I.I , 2) < radius)
                M[I] = 1
        end
    end
end
M
display(heatmap(M))
size(M)
ϕ = zeros(size(M) .+ (2,2))
ϕ[2:end-1 , 2:end-1] =  M
display(heatmap(ϕ))

discrete_G_weigted_neigbour_sum(19,12 , ϕ , __G_h , size(M ,1 ) , size(M ,2))

res = elyps_solver( zeros(size(M) .+ (2,2))  , ϕ , size(M , 1) , size(M , 2), 1000000001 , 1e-3 , 100000)
LinearAlgebra.norm(ϕ .- res)
display(heatmap(res))

W_prime(x) = -x * (1- x^2)

ξ , ψ = set_xi_and_psi(ϕ , size(M,1) , size(M,2), 1e-3 , W_prime)
display(heatmap(ψ))
μ = zeros(size(ψ))
μ[2:end-1,2:end-1] = W_prime.(ϕ[2:end-1,2:end-1])

display(heatmap(μ))




smoothed = SMOOTH(copy(ξ) , copy(ψ), copy(ϕ) , copy(μ) , 1e-3, 1e-3 , 1e-3 , size(M,1) , size(M,2) , 500 , true)
display(heatmap(smoothed[1]))
