include(pwd() * "/" * "utils.jl")

function L(solver::multi_solver,i,j , phi , mu)
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

function dL(solver::multi_solver , i , j)
    return [ (1/solver.dt) (1/solver.h^2*neighbours_in_domain(i,j,G,solver.len , solver.width));
             (-1*solver.epsilon^2/solver.h^2 * neighbours_in_domain(i,j,G,solver.len , solver.width) - 2) 1]
    end

function SMOOTH!(
    solver::multi_solver,
    iterations,
    adaptive
)
    for k = 1:iterations
        old_phase = copy(solver.phase)
        for I in CartesianIndices(solver.phase)[2:end-1, 2:end-1]
            i, j = I.I
            bordernumber = neighbours_in_domain(i, j, G, solver.len, solver.width)
            coefmatrix =
                [
                    (1/solver.dt) (bordernumber/solver.h^2);
                    (-1*(2+(solver.epsilon^2/solver.h^2)*bordernumber)) 1
                ]


            b =
                [
                    (
                        solver.xi[i, j]
                        +
                        discrete_G_weigted_neigbour_sum(
                            i, j, solver.potential, G, solver.len, solver.width
                        )
                        /
                        solver.h^2
                    ),
                    (
                        solver.psi[i, j]
                        -
                        (solver.epsilon^2 / solver.h^2)
                        *
                        discrete_G_weigted_neigbour_sum(
                            i, j, solver.phase, G, solver.len, solver.width
                        )
                    )
                ]

            res = coefmatrix \ b
            solver.phase[i, j] = res[1]
            solver.potential[i, j] = res[2]

        end

        if adaptive && LinearAlgebra.norm(old_phase - solver.phase) < 1e-8
            #println("SMOOTH terminated at $(k) succesfully")
            break
        end
    end
end

function v_cycle(grid::Array{multi_solver}, level)

    solver = grid[level]
    SMOOTH!(solver, 400, true)
    #println("Finished pre SMOOTHing")

    # extract (d,r) as array operations

    d = zeros(size(solver.phase))
    r = zeros(size(solver.phase))

    for I in CartesianIndices(solver.phase)[2:end-1, 2:end-1]
        d[I], r[I] = [solver.xi[I], solver.psi[I]] .- L(solver, I.I..., solver.phase[I], solver.potential[I])
    end

    # print(f"Max derivation d: {np.linalg.norm(d)}")
    # print(f"Max derivation r: {np.linalg.norm(r)}")
    restrict_solver!(grid[level], grid[level+1])
    solver = grid[level+1]
    solution = deepcopy(solver)

d_large = restrict(d, G)
r_large = restrict(r, G)

#println(" d $(norm(d_large))")
#println(" r $(norm(r_large))")

u_large = zeros(size(d_large))
v_large = zeros(size(d_large))

    #TODO short newton iteration for
    for i = 1:300
        for I in CartesianIndices(solver.phase)[2:end-1, 2:end-1]

            diffrence = L(solution, I.I..., solution.phase[I], solution.potential[I]) .- [d_large[I], r_large[I]] .- L(solver, I.I..., solver.phase[I], solver.potential[I])
            #diffrence = collect(L(solution, I.I...)) .- collect(L(solver, I.I...))
            #diffrence = [d_large[I] , r_large[I]]

            local ret = dL(solution, I.I...) \ diffrence

            u_large[I] = ret[1]
            v_large[I] = ret[2]
        end
        solution.phase .-= u_large
        solution.potential .-= v_large
    end
    #println("Finished  largegrid")

    u_large = solver.phase .- solution.phase
    v_large = solver.potential .- solution.potential

    solver = grid[level]

    solver.phase .+= prolong(u_large , G)
    solver.potential .+= prolong(v_large, G)
    # smooth again:
    SMOOTH!(solver, 800, true)
    #println("Finished post SMOOTHing")
end
