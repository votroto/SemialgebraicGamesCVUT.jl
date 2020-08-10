# -----------------------------------------------------------------------------
# Diagonally-dominant matrices [JuMP]
# -----------------------------------------------------------------------------

# Ahmadi, A., & Majumdar, A. (2017). Dsos and sdsos optimization: More
# tractable alternatives to sum of squares and semidefinite optimization.
# SIAM Journal on Applied Algebra and Geometry, 3. https://doi.org/
# 10.1137/18M118935X

# JuMP constraints for diagonal-dominance -------------------------------------

"""
Constraint for weak column diagonal dominance with nonnegative diagonal entries.
A matrix M is diagonally-dominant, or M ∈ dd, if:
        |M[i,i]| ≥ ∑_{i≂̸j} |M[i,j]|, ∀i.
However, for our problem, dd also imposes nonnegative diagonal entries.
"""
struct dd end
struct ddConstraint <: AbstractConstraint
        mat::AbstractArray
end

JuMP.build_constraint(::Function, m::AbstractArray, ::dd) = ddConstraint(m)

function JuMP.add_constraint(model::Model, c::ddConstraint, ::String = "")
        X = c.mat
        n = size(X, 1)
        W = X - Diagonal(X)
        Z = @variable(model, [1:n, 1:n], Symmetric)
        @constraint(model, W .>= -Z)
        @constraint(model, W .<= Z)
        @constraint(model, diag(X) .>= sum(Z, dims = 2))
end

# JuMP constraints for the DD family of cones ---------------------------------
# DD(U) := {M ∈ Symmetric | M = U'QU for some dd matrix Q}.
# One can think of DD(U) as the set of matrices that are dd after a change of
# coordinates via the matrix U.
struct DD
        u::AbstractArray
end
struct DDConstraint <: AbstractConstraint
        mat::AbstractArray
        u::AbstractArray
end

JuMP.build_constraint(::Function, m::AbstractArray, d::DD) =
        DDConstraint(m, d.u)

function JuMP.add_constraint(model::Model, c::DDConstraint, ::String = "")
        n = size(c.mat, 1)
        Q = @variable(model, [1:n, 1:n], Symmetric)
        @constraint(model, Q in dd())
        xs = JuMP.vectorize(c.mat - c.u' * Q * c.u, SymmetricMatrixShape(n))
	@constraint(model, xs in MOI.Zeros(length(xs)))
end