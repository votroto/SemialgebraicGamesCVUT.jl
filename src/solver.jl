using SumOfSquares: SOSModel
using MultivariateMoments
using MultivariatePolynomials

function min_order(p::AbstractPolynomial, Sy::AbstractSemialgebraicSet)
	# Finds the initial order of the Lasserre moment hierarchy; i.e what
	# moments are required to compute the expectations of p and Sy.

	round_up_even(i::Integer) = (i + 1) & ~(1)
	defs(s::AbstractBasicSemialgebraicSet) = inequalities(s)
	defs(s::AbstractAlgebraicSet) = equalities(s)

	round_up_even(maximum([maxdegree.(p, variables(Sy)); maxdegree.(defs(Sy))]))
end

function solve_game(p::AbstractPolynomial, Sx::AbstractSemialgebraicSet,
	Sy::AbstractSemialgebraicSet, optimizer; iteration::Integer = 0)
	# Tries to solve a zero-sum two-player polynomial game on semialgebraic
	# sets, using an iteration of the Lasserre moment hierarchy.

	# Returns the payoff and the optimal strategy of player X; given a
	# payoff function p, and semialgebraic strategy sets Sx and Sy.

	# Lasserre, Jean-Bernard. (2004). Global Optimization With Polynomials
	# And The Problem Of Moments. SIAM Journal on Optimization.
	# 11. 10.1137/S1052623400366802.
	# P. A. Parrilo (2006). Polynomial games and sum of squares optimization
	# In Proceedings of the 45th IEEE Conference on Decision and Control
	# (pp. 2855-2860).


	@assert iteration >= 0

	order = min_order(p, Sy) + iteration * 2
	monoms = reverse(monomials(variables(Sy), 0:order))

	m = SOSModel(optimizer)
	@variable  m α
	@variable  m μ[1:length(monoms)]
	@objective m Min α

	μs = measure(μ, monoms)

	@constraint(m, α - expect(μs, p) >= 0, domain = Sx)
	@constraint(m, μs in MomentSequence(), domain = Sy)
	@constraint(m, μ[1] == 1)

	set_silent(m)
	optimize!(m)

	mea = measure(value.(μ), monoms)
	mm = moment_matrix(mea, monomials(variables(Sy), 0:order ÷ 2))
	value(α), isnan(value(α)) ? nothing : extractatoms(mm, 1e-4)
end
