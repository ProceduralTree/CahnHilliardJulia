function set_xi_and_psi(phase , len , width , dt , W_prime)
    xi = zeros(len + 2 , width + 2)
    psi = zeros(len + 2 , width + 2)
    xi_init(x) = x / dt
    psi_init(x) = W_prime(x) - 2 * x
    xi[2:end-1, 2:end-1] = xi_init.(phase[2:end-1,2:end-1])
    psi[2:end-1, 2:end-1] = psi_init.(phase[2:end-1,2:end-1])

    (xi , psi)
end

function testdata(gridsize , blobs , radius , norm=2 )
rngpoints = rand(1:gridsize, 2, 10)
M = zeros(gridsize,gridsize) .- 1

for p  axes(rngpoints , 2)
    point = rngpoints[:, p]
    for I in eachindex(IndexCartesian(), M)
            if (LinearAlgebra.norm(point .- I.I  , 2) < radius)
                M[I] = 1
            end
    end
end
   return M
end
