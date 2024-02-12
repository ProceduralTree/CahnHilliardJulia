
include("/home/jon/Projects/julia_tst/utils.jl")
include("/home/jon/Projects/julia_tst/multisolver.jl")
using Plots
using LinearAlgebra
M = testdata(64, 10, 15, 2)

testgrd = testgrid(M, 2)
p0 = heatmap(testgrd[1].phase, title="Initial State");
solver = testgrd[1]
set_xi_and_psi!(solver)
SMOOTH!(solver, 200, true);
p1 = heatmap(solver.phase, title="After Pre Smoothing");


d = zeros(size(solver.phase))
r = zeros(size(solver.phase))

for I in CartesianIndices(solver.phase)[2:end-1, 2:end-1]
    d[I], r[I] = [solver.xi[I], solver.psi[I]] .- L(solver, I.I...)
end

p2 = heatmap(d, title="Phase Residuals");
level = 1

restrict_solver!(testgrd[level], testgrd[level+1])
solver =testgrd[level+1]
solution = deepcopy(solver)


d_large = zeros(size(solver.phase))
r_large = zeros(size(solver.phase))

d_large[2:end-1, 2:end-1] = restrict(d[2:end-1, 2:end-1])
r_large[2:end-1, 2:end-1] = restrict(r[2:end-1, 2:end-1])

println(" d $(norm(d_large))")
println(" r $(norm(r_large))")

for i = 1:10
    for I in CartesianIndices(solver.phase)[2:end-1, 2:end-1]

        #diffrence = L(solution, I.I...) .- [d_large[I], r_large[I]] .- L(solver, I.I...)
        #diffrence = collect(L(solution, I.I...)) .- collect(L(solver, I.I...))
        diffrence = [d_large[I] , r_large[I]]

        ret = dL(solution , I.I...) \ diffrence
        if I == CartesianIndex(32,32)  println("Diffrence: $(diffrence)") end

        solution.phase[I] += ret[1]
        solution.potential[I] += ret[2]
    end
end


u_large = solver.phase .- solution.phase
v_large = solver.potential .- solution.potential

# u_large , v_large = (solver.phase , solver.potential) .- (solution.phase , solution.potential)
p3 = heatmap(u_large, title="error $(norm(u_large)) ")
p = plot(p0, p1, p2,p3, layout=(2, 2));
savefig(p, "images/v_cycle.svg")
