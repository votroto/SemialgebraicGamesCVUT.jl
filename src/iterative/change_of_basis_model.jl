# WIP !!
const CI = MOI.ConstraintIndex
const VI = MOI.VariableIndex

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

`optimizer_constructor` requires optimizer constructors such as ECOS.Optimizer.
`dd_family`             is either the DD or SDD family of cones.
`iterations`            specifies the number of changes-of-basis performed.
`number_type`           sets the numeric type for calculation.

Ahmadi, A., & Majumdar, A. (2017). Dsos and sdsos optimization: More
tractable alternatives to sum of squares and semidefinite optimization.
SIAM Journal on Applied Algebra and Geometry, 3. https://doi.org/
10.1137/18M118935X
"""
function CoBModel(
        optimizer_constructor,
        dd_family::Type{<:Union{SDD,DD}};
        iterations::Unsigned = 0,
        number_type::Type = Float64,
)
        m = CoBModel{number_type}()
        m.ext[:optimizer_constructor] = optimizer_constructor
        m.ext[:iterations] = iterations
        m.ext[:dd_family] = dd_family
        m.ext[:silent] = false
        return m
end

function MOI.optimize!(data::CoBModel{T}) where {T}
        # Initial Us are identity matrices
        initial_us(sizes) = [Array{T}(I, n, n) for n in sizes]
        # Each following iteration decomposes the associated X matrices.
        cholesky_us(Xs) = [cholesky_like(value.(x)) for x in Xs]
        function set_side_dimension(backend, id)
                set = MOI.get(backend, MOI.ConstraintSet(), id)
                MOI.side_dimension(set)
        end
        # Retrieve JuMP matrix variable from a MOI upper triangle vector
        function matrix_variable(model, id)
                fun = MOI.get(backend(model), MOI.ConstraintFunction(), id)
                side = set_side_dimension(backend(model), id)
                square = square_a_triangle(fun, side)
                JuMP.VariableRef.(model, square)
        end
        # Convert PSD constraints to DD(U) or SDD(U).
        function relax_ddp(model, map, sdps, Us)
                relaxed_map = Dict{CI,CI}()
                Xs = Vector{Matrix{VariableRef}}(undef, length(sdps))
                for p in 1:length(sdps)
                        id = map[sdps[p]]
                        Xs[p] = matrix_variable(model, id)
                        dominance = @constraint(model, Xs[p] in family(Us[p]))
                        MOI.delete(backend(model), id)
                        relaxed_map[sdps[p]] = dominance.index
                end
                Xs, relaxed_map
        end

        psd_cis_type = MOI.ListOfConstraintIndices{
                MOI.VectorOfVariables,
                MOI.PositiveSemidefiniteConeTriangle,
        }()

        family = data.ext[:dd_family]
        sdp_ids = MOI.get(data, psd_cis_type)
        sdp_sizes = set_side_dimension.(data, sdp_ids)
        iterations = data.ext[:iterations]

        Us = initial_us(sdp_sizes)

        for i in 0:iterations
                # Convert SDP to DDP before attaching Optimizer
                model = Model()
                MOI.set(model, MOI.Silent(), data.ext[:silent])
                copy_map = MOI.copy_to(model, data)
                Xs, relaxed_map = relax_ddp(model, copy_map, sdp_ids, Us)
                JuMP.set_optimizer(model, data.ext[:optimizer_constructor])

                optimize!(model)

                # Save result references
                back = backend(model)
                data.ext[:inner_optimizer] = back
                data.ext[:psd_map] = relaxed_map
                data.ext[:model_to_optimizer_map] = copy_map

                if !available_result(back)
                        break
                end

                # Change basis
                Us = cholesky_us(Xs)
        end
end

const StatusTypes = Union{
        MOI.TerminationStatus,
        MOI.PrimalStatus,
        MOI.DualStatus,
        MOI.ObjectiveValue,
}

function MOI.get(model::CoBModel{T}, attr::StatusTypes) where {T}
        MOI.get(model.ext[:inner_optimizer], attr)
end

function MOI.get(model::CoBModel{T}, attr::MOI.VariablePrimal, vi::VI) where {T}
        idx = MOIU.map_indices(model.ext[:model_to_optimizer_map], vi)
        MOI.get(model.ext[:inner_optimizer], attr, idx)
end

function MOI.get(model::CoBModel{T}, attr::MOI.ConstraintDual, ci::CI) where {T}
        idx = MOIU.map_indices(model.ext[:model_to_optimizer_map], ci)
        MOI.get(model.ext[:inner_optimizer], attr, idx)
end

# Overload for constraints changed from PSD to DD or SDD.
function MOI.get(
        model::CoBModel{T},
        attr::MOI.ConstraintDual,
        ci::CI{F,MOI.PositiveSemidefiniteConeTriangle},
) where {T,F}
        idx = model.ext[:psd_map][ci]
        MOI.get(model.ext[:inner_optimizer], attr, idx)
end

MOI.supports(::CoBModel, ::MOI.Silent) = true
function MOI.set(optimizer::CoBModel, ::MOI.Silent, value::Bool)
        optimizer.ext[:silent] = value
end

Base.show(io::IO, ::CoBModel) = print(io, "An iterative change-of-basis model")
