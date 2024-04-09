function B(I::CartesianIndex,shape::Tuple)
    laplace = 1
    for (i,s) in zip(I.I , shape)
        if i == s || i == 1
            return laplace
        end
    end
    return 0
end

function L(solver::T,i,j , phi , mu) where T <: Union{adapted_multi_solver, adapted_relaxed_multi_solver , gradient_boundary_solver}
    xi = solver.phase[i, j] / solver.dt -
         (discrete_G_weigted_neigbour_sum(i, j, solver.potential, G, solver.len, solver.width)
          -
          neighbours_in_domain(i, j, G, solver.len, solver.width) * mu )/solver.h^2
    psi = solver.epsilon^2/solver.h^2 *
          (discrete_G_weigted_neigbour_sum(i, j, solver.phase, G, solver.len, solver.width)
           -
           neighbours_in_domain(i, j, G, solver.len, solver.width) * phi) - 2 * phi + mu
    return [xi, psi]
end

function set_xi_and_psi!(solver::T) where T <: Union{adapted_multi_solver, adapted_relaxed_multi_solver}
    xi_init(x) = x / solver.dt
    psi_init(x) = solver.W_prime(x) - 2 * x
    solver.xi[2:end-1, 2:end-1] = xi_init.(solver.phase[2:end-1,2:end-1])
    solver.psi[2:end-1, 2:end-1] = psi_init.(solver.phase[2:end-1,2:end-1]) + B.(CartesianIndices(solver.phase[2:end-1,2:end-1]) , Ref((solver.len , solver.width)) )
    return nothing
end

function dL(solver::T , i , j) where T <: Union{adapted_multi_solver, adapted_relaxed_multi_solver, gradient_boundary_solver}
    return [ (1/solver.dt) (1/solver.h^2*neighbours_in_domain(i,j,G,solver.len , solver.width));
             (-1*solver.epsilon^2/solver.h^2 * neighbours_in_domain(i,j,G,solver.len , solver.width) - 2) 1]
    end

function B_1(i,j,shape)
    B_x = 0.01
    B_y = 0
    return (G(i + 0.5,j, shape[1] , shape[2]) - G(i - 0.5,j,shape[1] , shape[2])) * B_x
         + (G(i,j + 0.5, shape[1] , shape[2]) - G(i,j - 0.5,shape[1] , shape[2])) * B_y
end

function div(f, I::CartesianIndex  , shape , h)
    (i,j) = I.I
    return (f(i+0.5 , j , shape) - f(i-0.5, j , shape) + f(i,j+0.5,shape) - f(i, j - 0.5, shape))/ h
    end
function set_xi_and_psi!(solver::T) where T <: Union{gradient_boundary_solver}
    xi_init(x) = x / solver.dt
    psi_init(x) = solver.W_prime(x) - 2 * x
    solver.xi[2:end-1, 2:end-1] = xi_init.(solver.phase[2:end-1,2:end-1])
    solver.psi[2:end-1, 2:end-1] = psi_init.(solver.phase[2:end-1,2:end-1]) + div.(B_1 ,CartesianIndices(solver.phase[2:end-1,2:end-1]) , Ref((solver.len , solver.width)) , Ref(solver.h) )
    return nothing
end
