# -----------------------------------------------------------------------------
# Scaled-diagonally-dominant matrices [JuMP]
# -----------------------------------------------------------------------------

# Ahmadi, A., & Majumdar, A. (2017). Dsos and sdsos optimization: More
# tractable alternatives to sum of squares and semidefinite optimization.
# SIAM Journal on Applied Algebra and Geometry, 3. https://doi.org/
# 10.1137/18M118935X

# JuMP constraints for scaled-diagonal-dominance ------------------------------
struct sdd end
struct sddConstraint <: AbstractConstraint
        mat::AbstractArray
end

JuMP.build_constraint(::Function, m::AbstractArray, ::sdd) = sddConstraint(m)

function JuMP.add_constraint(model::Model, c::sddConstraint, ::String = "")
        X = Matrix{GenericAffExpr{Float64,VariableRef}}(c.mat)
        n = size(X, 1)
        for j in 1:n, i in 1:j-1
                r = @variable(model, [1:3])
                g = [r[1], r[3], √2 * r[2]]
                @constraint(model, g in MOI.RotatedSecondOrderCone(3))
                X[i, i] -= r[1]
                X[i, j] -= r[2]
                X[j, i] -= r[2]
                X[j, j] -= r[3]
        end
        @constraint(model, X .== zeros(size(X)))
end

# JuMP constraints for the SDD family of cones --------------------------------
# SDD(U) := {M ∈ Symmetric | M = U'QU for some sdd matrix Q}.
# One can think of SDD(U) as the set of matrices that are sdd after a change of
# coordinates via the matrix U.
struct SDD
        u::AbstractArray
end
struct SDDConstraint <: AbstractConstraint
        mat::AbstractArray
        u::AbstractArray
end

JuMP.build_constraint(::Function, m::AbstractArray, d::SDD) =
        SDDConstraint(m, d.u)

function JuMP.add_constraint(model::Model, c::SDDConstraint, ::String = "")
        n = size(c.mat, 1)
        Q = @variable(model, [1:n, 1:n], Symmetric)
        @constraint(model, Q in sdd())
        xs = JuMP.vectorize(c.mat - c.u' * Q * c.u, SymmetricMatrixShape(n))
        @constraint(model, xs in MOI.Zeros(length(xs)))
end
