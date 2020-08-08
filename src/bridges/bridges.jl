using MathOptInterface
const MOI = MathOptInterface
const MOIB = MathOptInterface.Bridges
const MOIU = MathOptInterface.Utilities

function square_a_triangle(f::MOI.AbstractVectorFunction, dimension::Integer)
        x = MOIU.eachscalar(f)
        S = MOIU.scalar_type(typeof(f))
        matrix = Matrix{S}(undef, dimension, dimension)
        k = 0
        for j in 1:dimension, i in 1:j
                k += 1
                matrix[j, i] = matrix[i, j] = x[k]
        end
        matrix
end

include("sddp_bridge.jl")
include("ddp_bridge.jl")
include("lazy_relax.jl")
