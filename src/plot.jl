include(pwd() * "/src/" * "solvers.jl")
include(pwd() * "/src/" * "adapted_solvers.jl")
include(pwd() * "/src/" * "utils.jl")
include(pwd() * "/src/" * "multisolver.jl")
include(pwd() * "/src/" * "testgrids.jl")
using Plots
using LaTeXStrings
using LinearAlgebra
using Printf
M = testdata(64, 16, 12 , 2)

testgrd = testgrid(multi_solver,M, 2)
test_solver = testgrd[1]

p0 = heatmap(testgrd[1].phase, title="Initial State");
solver = testgrd[1]
set_xi_and_psi!(solver)
SMOOTH!(solver, 400, true);
p1 = heatmap(solver.phase, title="After Pre Smoothing");


d = zeros(size(solver.phase))
r = zeros(size(solver.phase))

for I in CartesianIndices(solver.phase)[2:end-1, 2:end-1]
    d[I], r[I] = [solver.xi[I], solver.psi[I]] .- L(solver, I.I..., solver.phase[I] , solver.potential[I])
end

p2 = heatmap(d, title="Phase Residuals");
level = 1

restrict_solver!(testgrd[level], testgrd[level+1])
solver =testgrd[level+1]
solution = deepcopy(solver)



d_large = restrict(d, G)
r_large = restrict(r, G)

println(" d $(norm(d_large))")
println(" r $(norm(r_large))")

u_large = zeros(size(d_large))
v_large = zeros(size(d_large))



for i = 1:300
    for I in CartesianIndices(solver.phase)[2:end-1, 2:end-1]


        diffrence = L(solution, I.I..., solution.phase[I], solution.potential[I]) .- [d_large[I], r_large[I]] .- L(solver, I.I... , solver.phase[I] , solver.potential[I])
        #diffrence = collect(L(solution, I.I...)) .- collect(L(solver, I.I...))
        #diffrence = [d_large[I] , r_large[I]]

        local ret = dL(solution , I.I...) \ diffrence
        #if I == CartesianIndex(2,2)  println("Diffrence: $(diffrence) , Ret: $(ret)") end

        u_large[I] = ret[1]
        v_large[I] = ret[2]
    end
    solution.phase .-= u_large
    solution.potential .-= v_large
end


p3 = heatmap(u_large, title=@sprintf "Change: %.1e" norm(u_large))
p = plot(p0, p1, p2,p3, layout=(2, 2));
savefig(p, "images/v_cycle.svg")

include(pwd() * "/src/" * "solvers.jl")
include(pwd() * "/src/" * "adapted_solvers.jl")
include(pwd() * "/src/" * "utils.jl")
include(pwd() * "/src/" * "multisolver.jl")
include(pwd() * "/src/" * "testgrids.jl")
using Plots
using LaTeXStrings
using LinearAlgebra
using Printf
M = testdata(64, 16, 12 , 2)

testgrd = testgrid(multi_solver,M, 2)
test_solver = testgrd[1]
set_xi_and_psi!(solver)

pbar = ProgressBar(total = 1000)

anim = @animate for i in 1:100
    for j in 1:10
        v_cycle!(testgrd, 1)
        update(pbar)
        end
    set_xi_and_psi!(testgrd[1])
    heatmap(testgrd[1].phase , clim =(-1,1) , framestyle=:none )
end
gif(anim , "images/iteration.gif" , fps = 10)
