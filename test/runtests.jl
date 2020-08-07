using SemialgebraicGames
using Test

tests = ["system_univariate.jl", "system_multivariate.jl"]

for t in tests
    include(t)
end
