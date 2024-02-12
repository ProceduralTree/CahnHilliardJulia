using RestrictProlong


"""
    restrict!(smallgrid_solver::multi_solver , largegrid_solver::multi_solver)::multi_solver

------------
Requires
----------
smallgrid solver and largegid solvers to be multiple of 2 from each other bar padding eg. (66x66)->(34x34)


TBW
"""
function restrict_solver!(smallgrid_solver::multi_solver, largegrid_solver::multi_solver)
    largegrid_solver.phase[2:end-1, 2:end-1] = restrict(smallgrid_solver.phase[2:end-1, 2:end-1])
    largegrid_solver.potential[2:end-1, 2:end-1] = restrict(smallgrid_solver.potential[2:end-1, 2:end-1])
    largegrid_solver.xi[2:end-1, 2:end-1] = restrict(smallgrid_solver.xi[2:end-1, 2:end-1])
    largegrid_solver.psi[2:end-1, 2:end-1] = restrict(smallgrid_solver.psi[2:end-1, 2:end-1])
    return nothing
end
function prolong!(smallgrid::multi_solver, largegrid::multi_solver)

end
