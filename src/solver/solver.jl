using SumOfSquares
using MultivariateMoments
using MultivariatePolynomials
using PolyJuMP

"""
Finds the initial order of the hierarchy; i.e what degrees are required
to compute: Eμ[p]; M[gμ] ∀g ∈ Sy; and positivstellensatz of Eμ[p] on Sx.
"""
function min_order(
        p::AbstractPolynomial,
        Sx::AbstractSemialgebraicSet,
        Sy::AbstractSemialgebraicSet,
)
        round_up_even(i::Integer) = (i + 1) & ~(1)
        deg(s::AbstractBasicSemialgebraicSet) = maxdegree.(inequalities(s))
        deg(s::AbstractAlgebraicSet) = maxdegree.(equalities(s))
        deg(p::AbstractPolynomial) = maxdegree(p)

        y = variables(Sy)
        x = variables(Sx)
        px = subs(p, y => ones(length(y)))
        py = subs(p, x => ones(length(x)))

        round_up_even(maximum([deg(px); deg(py); deg(Sy); deg(Sx)]))
end

"""
Tries to solve a zero-sum two-player polynomial game on semialgebraic sets,
using an iteration of the moment--sos hierarchy.

The inputs are: a payoff function p, semialgebraic strategy sets Sx and Sy, and
an "iteration", where "minimal order + 2*iteration" defines the actual order.

Returns: the payoff and the optimal strategies of both players.

Lasserre, Jean-Bernard. (2004). Global Optimization With Polynomials And The
Problem Of Moments. SIAM Journal on Optimization. 11. 10.1137/S1052623400366802.
P. A. Parrilo (2006). Polynomial games and sum of squares optimization In
Proceedings of the 45th IEEE Conference on Decision and Control (pp. 2855-2860).
"""
function solve_game(
        p::AbstractPolynomial,
        Sx::AbstractSemialgebraicSet,
        Sy::AbstractSemialgebraicSet,
        optimizer_factory;
        iteration::Unsigned = 0x0,
)
        order = min_order(p, Sx, Sy) + iteration * 2

        m = SOSModel(optimizer_factory)

        @variable m α
        @objective m Min α
        @variable m μ Moments(order, Sy)
        @constraint(m, q, linearize(μ, p) <= α, domain = Sx, maxdegree = order)

        set_silent(m)
        optimize!(m)

        if !available_result(m)
                return JuMP.termination_status(m), nothing, nothing
        end

        xmat = moment_matrix(q)
        ymat = moment_matrix(μ)
        value(α), extractatoms(xmat, 1e-4), extractatoms(ymat, 1e-4)
end
