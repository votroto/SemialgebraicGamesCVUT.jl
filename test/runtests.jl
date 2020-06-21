using SemialgebraicGames
using Test

tests = ["test_solver.jl"]

for t in tests
    include(t)
end
