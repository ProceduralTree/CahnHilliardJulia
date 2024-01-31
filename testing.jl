


include("elypssolver.jl")

using Plots
using  LinearAlgebra
SIZE = 2^9
radius = 10
rngpoints = rand(1:SIZE, 2, 10)
M = zeros(SIZE,SIZE) .- 1



for p = 1:size(rngpoints, 2)
    point = rngpoints[:, p]
    for i = 1:size(M , 1)
        for j = 1:size(M , 2)
            if (LinearAlgebra.norm(point .- [i,j] , 2) < radius)
                M[i,j] = 1
           end
        end
    end
end
M
display(heatmap(M))
size(M)
phase = zeros(size(M) .+ (2,2))
phase[2:end-1 , 2:end-1] =  M
display(heatmap(phase))

discrete_G_weigted_neigbour_sum(19,12 , phase , __G_h , size(M ,1 ) , size(M ,2))

res = elyps_solver( zeros(size(M) .+ (2,2))  , phase , size(M , 1) , size(M , 2), 1000000001 , 1e-3 , 100000)
LinearAlgebra.norm(phase .- res)
display(heatmap(res))
