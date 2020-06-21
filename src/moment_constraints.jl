using MultivariatePolynomials: AbstractPolynomialLike
using MultivariateMoments: AbstractMeasure
using LinearAlgebra: dot
using SemialgebraicSets
using JuMP

# Pirates an AlgebraicSet from SemialgebraicSets as a workaround
# for a design problem.
import SemialgebraicSets.inequalities
function SemialgebraicSets.inequalities(::AlgebraicSet)
	[]
end

struct MomentSequence end
struct MomentSequenceConstraint <: AbstractConstraint
	measure::AbstractMeasure
	domain::AbstractSemialgebraicSet
end

function expect(m::AbstractMeasure, p::AbstractPolynomialLike)
	# Comptes ∫p dμ = ∑(c * μ[i]),
	# for a given p = ∑(c * x^i).

	vars_p = variables(p)
	vars_m = variables(m)
	vars_d = setdiff(vars_p, vars_m)
	vec_coeffs = subs.(terms(p), vars_m => ones(size(vars_m)))
	vec_monoms = subs.(monomials(p), vars_d => ones(size(vars_d)))
	vec_coeffs' * (ms->dot(m, ms)).(vec_monoms)
end

JuMP.build_constraint(::Function, m::AbstractMeasure, ::MomentSequence; domain) =
	MomentSequenceConstraint(m, domain)

function JuMP.add_constraint(model::Model, c::MomentSequenceConstraint, ::String = "")
	# Builds Lasserre-style moment constraint of order 't'.
	# Mt(μ) >= 0; Mt-g(gμ) >= 0 (∀g); Mt-h(hμ) == 0 (∀h)

	# Lasserre, Jean-Bernard. (2004). Global Optimization With Polynomials
	# And The Problem Of Moments. SIAM Journal on Optimization.
	# 11. 10.1137/S1052623400366802.

	function moment_matrix(m::AbstractMeasure, t::Integer)
		ys = reverse(monomials(variables(m), 0:t ÷ 2))
		(x->expect(m, x)).(ys * ys')
	end
	function moment_matrix(m::AbstractMeasure, t::Integer, g::AbstractPolynomialLike)
		round_up_even(i::Integer) = (i + 1) & ~(1)
		tg = round_up_even(maxdegree(g))
		ys = reverse(monomials(variables(m), 0:t ÷ 2 - tg ÷ 2))
		(x->expect(m, x)).(ys * ys' * g)
	end

	t = maxdegree(monomials(c.measure))

	Mt = moment_matrix(c.measure, t)
	Mi = (p->moment_matrix(c.measure, t, p)).(inequalities(c.domain))
	Me = (p->moment_matrix(c.measure, t, p)).(equalities(c.domain))

	@constraint(model, Mt in PSDCone())
	@constraint(model, [i in Mi], i in PSDCone())
	@constraint(model, [i in Me], i .== zeros(size(i)))
end
