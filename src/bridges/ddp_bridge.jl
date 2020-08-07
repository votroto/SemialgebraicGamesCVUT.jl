"""
Converts a Semidefinite Program (SDP) into a Diagonally-Dominant Program (DDP).
Positive-semidefiniteness constraints are relaxed into Diagonal-dominance
constraints.
"""
struct DDPBridge{T,F,G} <: MOIB.Constraint.AbstractBridge end

function MOIB.Constraint.bridge_constraint(
        ::Type{DDPBridge{T,F,G}},
        m::MOI.ModelLike,
        f::MOI.AbstractVectorFunction,
        s::MOI.PositiveSemidefiniteConeTriangle,
) where {T,F,G}
        n = s.side_dimension
        z = MOI.add_variables(m, n)
        v = MOI.SingleVariable.(z)
        matrix = square_a_triangle(f, n)

        for i in 1:n
                zs = MOIU.operate(vcat, T, v[i], matrix[i, :]...)
                dominance = T(2) * matrix[i, i] - v[i]
                MOI.add_constraint(m, zs, MOI.NormOneCone(n + 1))
                MOI.add_constraint(m, dominance, MOI.GreaterThan(zero(T)))
        end
        return DDPBridge{T,F,G}()
end

function MOI.supports_constraint(
        ::Type{<:DDPBridge},
        ::Type{<:MOI.AbstractVectorFunction},
        ::Type{<:MOI.PositiveSemidefiniteConeTriangle},
)
        return true
end
function MOIB.added_constrained_variable_types(::Type{<:DDPBridge})
        return Tuple{DataType}[]
end
function MOIB.added_constraint_types(::Type{<:DDPBridge{T,F,G}}) where {T,F,G}
        return [(F, MOI.GreaterThan{T}); (G, MOI.NormOneCone)]
end
function MOIB.Constraint.concrete_bridge_type(
        ::Type{<:DDPBridge{T}},
        G::Type{<:MOI.AbstractVectorFunction},
        ::Type{MOI.PositiveSemidefiniteConeTriangle},
) where {T}
        S = MOIU.scalar_type(G)
        F = MOIU.promote_operation(-, T, S, MOI.SingleVariable)
        return DDPBridge{T,F,G}
end
