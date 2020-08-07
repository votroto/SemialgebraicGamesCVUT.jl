using MultivariatePolynomials
using DynamicPolynomials
using SemialgebraicSets
using JuMP

import MultivariateMoments.moment_matrix

struct Moments
        maxdegree::Integer
        domain::AbstractSemialgebraicSet
end
struct MomentsVar{T}
        maxdegree::Integer
        pvars::Vector{PolyVar{T}}
        mvars::AbstractArray{VariableRef,1}
        model
end
function JuMP.build_variable(_error, info, s::Moments; extra_kw_args...)
        return s
end

function moment_matrix(m::MomentsVar)
	ys = reverse(monomials(m.pvars, 0:m.maxdegree÷2))
	MomentMatrix(value.(linearize.(m, ys * ys')), ys)
end
function moment_matrix(m::MomentsVar, t::Integer)
	ys = reverse(monomials(m.pvars, 0:t÷2))
	linearize.(m, ys * ys')
end
function moment_matrix(m::MomentsVar, t::Integer, g::Polynomial)
	round_up_even(i::Integer) = (i + 1) & ~(1)
	tg = round_up_even(maxdegree(g))
	ys = reverse(monomials(m.pvars, 0:t÷2-tg÷2))
	linearize.(m, ys * ys' * g)
end

function JuMP.add_variable(model, s::Moments, name::String)
        pvars = variables(s.domain)
        monos = monomials(pvars, 0:s.maxdegree)
        mvars = @variable(model, [monos], base_name = name)
        mseq = MomentsVar(s.maxdegree, pvars, mvars, model)

        Mt = moment_matrix(mseq, s.maxdegree)
        Mi = moment_matrix.(mseq, s.maxdegree, inequalities(s.domain))
        Me = moment_matrix.(mseq, s.maxdegree, equalities(s.domain))

        @constraint(model, mvars[end] == 1)
        @constraint(model, Mt in PSDCone())
        @constraint(model, [i in Mi], i in PSDCone())
        @constraint(model, [i in Me], i .== zero(i))

        return mseq
end

linearize(m::MomentsVar, t::Term) = coefficient(t) * linearize(m, monomial(t))
linearize(m::MomentsVar, t::Polynomial) = sum(linearize.(m, terms(t)))
function linearize(m::MomentsVar, mono::Monomial)
	# workaround DynamicPolynomials' type instability

        if m.pvars == variables(mono)
                m.mvars[mono]
        else
                x = m.pvars
                ty = subs(mono, x => ones(length(x)))
                y = variables(ty)
                tx = subs(mono, y => ones(length(y)))

                ty * m.mvars[tx]
        end
end

Base.broadcastable(x::MomentsVar) = Ref(x)
