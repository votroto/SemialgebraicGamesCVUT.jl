using DynamicPolynomials
using MultivariateMoments
using SemialgebraicSets
using ECOS

# -----------------------------------------------------------------------------
# Trivial tests
# -----------------------------------------------------------------------------

@testset "Guessing game --------- " begin
	@polyvar x y
	p  = (x - y)^2
	Sx = @set 1 - x^2 >= 0
	Sy = @set 1 - y^2 >= 0
	opt = lazy_relax(ECOS.Optimizer, SDDPBridge)

	expected_measure_x = AtomicMeasure(variables(Sx), [
		WeightedDiracMeasure([1], 0.5),
		WeightedDiracMeasure([-1], 0.5)])
	expected_measure_y = AtomicMeasure(variables(Sy), [
		WeightedDiracMeasure([0], 1)])
	expected_value = 1

	pvalue, pstrat_x, pstrat_y = solve_game(+p, Sx, Sy, opt, iteration=0x1)
	dvalue, dstrat_y, dstrat_x = solve_game(-p, Sy, Sx, opt, iteration=0x1)

	eps = 1e-3
	@test isapprox(expected_value, pvalue, atol=eps)
	@test isapprox(pvalue, -dvalue, atol=eps)
	@test isapprox(expected_measure_x, dstrat_x, atol=eps)
	@test isapprox(expected_measure_y, pstrat_y, atol=eps)
end