# -----------------------------------------------------------------------------
# Scaled-Diagonally-dominant bridge [MOI]
# -----------------------------------------------------------------------------

"""
Converts a Semidefinite Program (SDP) into a Scaled=Diagonally-Dominant Program
(SDDP). Positive-semidefiniteness constraints are relaxed into Scaled-diagonal-
dominance constraints.
"""
struct SDDPBridge{T,F,G} <: MOIB.Constraint.AbstractBridge
        side_dimension
        constraint
end

function MOIB.Constraint.bridge_constraint(
        ::Type{SDDPBridge{T,F,G}},
        model::MOI.ModelLike,
        f::MOI.AbstractVectorFunction,
        s::MOI.PositiveSemidefiniteConeTriangle,
) where {T,F,G}
        n = s.side_dimension
        matrix = square_a_triangle(f, n)

        for j in 1:n, i in 1:j-1
                zs = MOI.add_variables(model, 3)
                vs = MOI.SingleVariable.(zs)
                MOIU.operate!(-, T, matrix[i, i], vs[1])
                MOIU.operate!(-, T, matrix[i, j], vs[2])
                MOIU.operate!(-, T, matrix[j, j], vs[3])

                g = MOIU.operate(vcat, T, vs[1], vs[3], T(√2) * vs[2])
                l = MOI.add_constraint(model, g, MOI.RotatedSecondOrderCone(3))
        end

        # Handle accidental PSD constraints on a single variable
        e = undef
        if n == 1
                g = matrix[1]
                e = MOI.add_constraint(model, g, MOI.GreaterThan(zero(T)))
        else
                g = MOIU.operate(vcat, T, matrix...)
                e = MOI.add_constraint(model, g, MOI.Zeros(n * n))
        end

        return SDDPBridge{T,F,G}(n, e)
end

function MOI.supports_constraint(
        ::Type{<:SDDPBridge},
        ::Type{<:MOI.AbstractVectorFunction},
        ::Type{<:MOI.PositiveSemidefiniteConeTriangle},
)
        return true
end

function MOIB.added_constrained_variable_types(::Type{<:SDDPBridge})
        return Tuple{DataType}[]
end
function MOIB.added_constraint_types(::Type{<:SDDPBridge{T,F,G}}) where {T,F,G}
        return [
                (G, MOI.RotatedSecondOrderCone)
                (G, MOI.Zeros)
                (F, MOI.GreaterThan{T})
        ]
end
function MOIB.Constraint.concrete_bridge_type(
        ::Type{<:SDDPBridge{T}},
        G::Type{<:MOI.AbstractVectorFunction},
        ::Type{MOI.PositiveSemidefiniteConeTriangle},
) where {T}
        S = MOIU.scalar_type(G)
        F = MOIU.promote_operation(-, T, S, MOI.SingleVariable)
        return SDDPBridge{T,F,G}
end

function MOI.get(ml::MOI.ModelLike, cd::MOI.ConstraintDual, bridge::SDDPBridge)
        dual = MOI.get(ml, cd, bridge.constraint)
        side = bridge.side_dimension
        mat = reshape(dual, side, side)
        [(i == j) ? mat[i, j] : mat[i, j] / 2 for j in 1:side for i in 1:j]
end
