
function set_xi_and_psi(phase , len , width , dt , W_prime)
    xi = zeros(len + 2 , width + 2)
    psi = zeros(len + 2 , width + 2)
    xi_init(x) = x / dt
    psi_init(x) = W_prime(x) - 2 * x
    xi[2:end-1, 2:end-1] = xi_init.(phase[2:end-1,2:end-1])
    psi[2:end-1, 2:end-1] = psi_init.(phase[2:end-1,2:end-1])

    (xi , psi)
end
