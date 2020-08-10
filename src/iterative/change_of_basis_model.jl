# WIP !!

MOI.Utilities.@model(
        CoBModel,
        (),
        (MOI.EqualTo, MOI.GreaterThan, MOI.LessThan),
        (MOI.Nonnegatives, MOI.PositiveSemidefiniteConeTriangle),
        (),
        (),
        (MOI.ScalarAffineFunction,),
        (MOI.VectorOfVariables,),
        (),
        true,
)

#function MOI.supports_constraint(
#        ::CoBModel{T},
#        ::Type{MOI.SingleVariable},
#        ::Type{<:MOI.Utilities.SUPPORTED_VARIABLE_SCALAR_SETS{T}},
#) where {T}
#
#        return false
#end

function MOI.supports(
        ::CoBModel{T},
        ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}},
) where {T}
        return true
end

function MOI.supports_constraint(
        ::CoBModel,
        ::Type{MOI.VectorOfVariables},
        ::Type{<:Union{MOI.Nonnegatives,MOI.PositiveSemidefiniteConeTriangle}},
)
        return true
end
function MOI.supports_constraint(
        ::CoBModel{T},
        ::Type{MOI.ScalarAffineFunction{T}},
        ::Type{<:Union{MOI.EqualTo{T},MOI.GreaterThan{T},MOI.LessThan{T}}},
) where {T}
        return true
end

"""
A change-of-basis iterative model.

Implements an improvement on DDP and SDDP optimization using the DD and SDD
families of cones.

Ahmadi, A., & Majumdar, A. (2017). Dsos and sdsos optimization: More
tractable alternatives to sum of squares and semidefinite optimization.
SIAM Journal on Applied Algebra and Geometry, 3. https://doi.org/
10.1137/18M118935X
"""
function CoBModel(
        optimizer_constructor,
        diagonal_dominance_family::Type{<:Union{SDD,DD}};
        iterations::Unsigned = 0,
        number_type::Type = Float64,
)
        m = CoBModel{number_type}()
        m.ext[:diagonal_dominance_family] = diagonal_dominance_family
        m.ext[:optimizer_constructor] = optimizer_constructor
        m.ext[:iterations] = iterations
        return m
end

function MOI.optimize!(data::CoBModel{T}) where {T}
        psd_vec = MOI.ListOfConstraintIndices{
                MOI.VectorOfVariables,
                MOI.PositiveSemidefiniteConeTriangle,
        }()

        ci = MOI.get(data, psd_vec)
        num_psd = length(ci)
        Us = Vector{Matrix{T}}(undef, num_psd)

        js = Vector{Matrix{VariableRef}}(undef, num_psd)
        for i in 0:data.ext[:iterations]
                m = Model()
                intermediate_map = MOI.copy_to(m, data)

                data.ext[:psd_map] =
                        Dict{MOI.ConstraintIndex,MOI.ConstraintIndex}()

                back = backend(m)

                ci = MOI.get(back, psd_vec)
                for p in 1:num_psd
                        fun = MOI.get(back, MOI.ConstraintFunction(), ci[p])
                        side = MOI.side_dimension(MOI.get(
                                back,
                                MOI.ConstraintSet(),
                                ci[p],
                        ))
                        U = isassigned(Us, p) ? Us[p] : Array{T}(I, side, side)
                        sqr = square_a_triangle(fun, side)
                        js[p] = JuMP.VariableRef.(m, sqr)
                        xxx = @constraint(m, js[p] in data.ext[:diagonal_dominance_family](U))
                        MOI.delete(back, ci[p])
                        data.ext[:psd_map][ci[p]] = xxx.index
                end
                JuMP.set_optimizer(m, data.ext[:optimizer_constructor])

                if data.ext[:silent]
                        JuMP.set_silent(m)
                end
                optimize!(m)

                data.ext[:inner_optimizer] = back
                data.ext[:model_to_optimizer_map] = intermediate_map
                for p in 1:num_psd
                        Us[p] = cholesky_like(value.(js[p]))
                end
        end
end

function MOI.get(model::CoBModel{T}, attr::MOI.TerminationStatus) where {T}
        MOI.get(model.ext[:inner_optimizer], attr)
end

function MOI.get(model::CoBModel{T}, attr::MOI.ObjectiveValue) where {T}
        MOI.get(model.ext[:inner_optimizer], attr)
end

function MOI.get(
        model::CoBModel{T},
        attr::MOI.VariablePrimal,
        vi::MOI.VariableIndex,
) where {T}
        MOI.get(
                model.ext[:inner_optimizer],
                attr,
                MOIU.map_indices(model.ext[:model_to_optimizer_map], vi),
        )
end

function MOI.get(
        model::CoBModel{T},
        attr::MOI.ConstraintDual,
        vi::MOI.ConstraintIndex,
) where {T}

        if haskey(model.ext[:psd_map], model.ext[:model_to_optimizer_map][vi])
                return MOI.get(
                        model.ext[:inner_optimizer],
                        attr,
                        model.ext[:psd_map][model.ext[:model_to_optimizer_map][vi]],
                )
        else
                return MOI.get(
                        model.ext[:inner_optimizer],
                        attr,
                        MOIU.map_indices(
                                model.ext[:model_to_optimizer_map],
                                vi,
                        ),
                )
        end
end

MOI.supports(::CoBModel, ::MOI.Silent) = true
function MOI.set(optimizer::CoBModel, ::MOI.Silent, value::Bool)
        optimizer.ext[:silent] = value
end
