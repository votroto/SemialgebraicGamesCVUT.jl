using MathOptInterface
const MOI = MathOptInterface

"""MOI utility to create a square matrix from an upper-triangle vector"""
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