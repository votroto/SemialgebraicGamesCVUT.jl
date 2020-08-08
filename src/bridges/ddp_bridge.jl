"""
Converts a Semidefinite Program (SDP) into a Diagonally-Dominant Program (DDP).
Positive-semidefiniteness constraints are relaxed into Diagonal-dominance
constraints.
"""
struct DDPBridge{T,F,G} <: MOIB.Constraint.AbstractBridge
        dominance::Vector{MOI.ConstraintIndex{F,MOI.EqualTo{T}}}
end

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

        ds = Vector{MOI.ConstraintIndex{F,MOI.EqualTo{T}}}(undef, n)
        for i in 1:n
                zs = MOIU.operate(vcat, T, v[i], matrix[i, :]...)
                dominance = T(2) * matrix[i, i] - v[i]
                MOI.add_constraint(m, zs, MOI.NormOneCone(n + 1))
                ds[i] = MOI.add_constraint(m, dominance, MOI.EqualTo(zero(T)))
        end
        return DDPBridge{T,F,G}(ds)
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
        return [(F, MOI.EqualTo{T}); (G, MOI.NormOneCone)]
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

function MOI.get(
        model::MOI.ModelLike,
        attr::MOI.ConstraintDual,
        bridge::DDPBridge{T},
) where {T}
        diag_dual = MOI.get(model, attr, bridge.dominance)
        side = length(diag_dual)
        n = MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(side))

        dd_dual = Vector{T}(undef, n)
        k = 0
        for j in 1:side
                for i in 1:(j-1)
                        k += 1
                        dd_dual[k] = -(diag_dual[i] + diag_dual[j])
                end
                k += 1
                dd_dual[k] = diag_dual[j] * 2
        end
        return dd_dual
end
