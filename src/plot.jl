include(pwd() * "/src/solvers.jl")
include(pwd() * "/src/adapted_solvers.jl")
include(pwd() * "/src/utils.jl")
include(pwd() * "/src/multisolver.jl")
include(pwd() * "/src/multi_relaxed.jl")
include(pwd() * "/src/testgrids.jl")
include(pwd() * "/src/elypssolver.jl")
using Plots
using LaTeXStrings
using LinearAlgebra
using Printf
using ProgressBars
default(fontfamily="computer modern")
SIZE = 64
M = testdata(SIZE, SIZE ÷ 5, SIZE /5 , 2)
testgrd = testgrid(multi_solver,M, 2)
test_solver = testgrd[1]

p0 = heatmap(testgrd[1].phase, title="Initial State");
s = testgrd[1]
set_xi_and_psi!(s)
SMOOTH!(s, 400, true);
p1 = heatmap(s.phase, title="After Pre Smoothing");


d = zeros(size(s.phase))
r = zeros(size(s.phase))

for I in CartesianIndices(s.phase)[2:end-1, 2:end-1]
    d[I], r[I] = [s.xi[I], s.psi[I]] .- L(s, I.I..., s.phase[I] , s.potential[I])
end

p2 = heatmap(d, title="Phase Residuals");
level = 1

restrict_solver!(testgrd[level], testgrd[level+1])
s =testgrd[level+1]
solution = deepcopy(s)



d_large = restrict(d, G)
r_large = restrict(r, G)

println(" d $(norm(d_large))")
println(" r $(norm(r_large))")

u_large = zeros(size(d_large))
v_large = zeros(size(d_large))



for i = 1:300
    for I in CartesianIndices(s.phase)[2:end-1, 2:end-1]


        diffrence = L(solution, I.I..., solution.phase[I], solution.potential[I]) .- [d_large[I], r_large[I]] .- L(s, I.I... , s.phase[I] , s.potential[I])
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

include(pwd() * "/src/solvers.jl")
include(pwd() * "/src/adapted_solvers.jl")
include(pwd() * "/src/utils.jl")
include(pwd() * "/src/multisolver.jl")
include(pwd() * "/src/multi_relaxed.jl")
include(pwd() * "/src/testgrids.jl")
include(pwd() * "/src/elypssolver.jl")
using Plots
using LaTeXStrings
using LinearAlgebra
using Printf
using ProgressBars
default(fontfamily="computer modern")
SIZE = 64
M = testdata(SIZE, SIZE ÷ 5, SIZE /5 , 2)
using JLD2
using DataFrames
results = jldopen("experiments/iteration.jld2")["result"]
anim = @animate for res in eachrow(results)
    heatmap(res.solver.phase , xlims = (2,size(res.solver.phase , 1)-1) , ylim=(2,size(res.solver.phase , 1)-1) , aspectratio=:equal)
end
gif(anim , "images/iteration.gif" , fps = 10)
