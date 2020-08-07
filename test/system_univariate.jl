using DynamicPolynomials
using MultivariateMoments
using SemialgebraicSets
using CSDP

# -----------------------------------------------------------------------------
# Trivial tests
# -----------------------------------------------------------------------------

@testset "Guessing game --------- " begin
	@polyvar x y
	p  = (x - y)^2
	Sx = @set 1 - x^2 >= 0
	Sy = @set 1 - y^2 >= 0
	opt = CSDP.Optimizer

	expected_measure_x = AtomicMeasure(variables(Sx), [
		WeightedDiracMeasure([1], 0.5),
		WeightedDiracMeasure([-1], 0.5)])
	expected_measure_y = AtomicMeasure(variables(Sy), [
		WeightedDiracMeasure([0], 1)])
	expected_value = 1

	value, measure_x, measure_y = solve_game(p, Sx, Sy, opt, iteration=0x1)

	eps = 1e-3
	@test isapprox(expected_value, value, atol=eps)
	@test isapprox(expected_measure_x, measure_x, atol=eps)
	@test isapprox(expected_measure_y, measure_y, atol=eps)
end

# -----------------------------------------------------------------------------
# P. A. Parrilo (2006). Polynomial games and sum of squares optimization.
# In Proceedings of the 45th IEEE Conference on Decision and Control
# (pp. 2855-2860).
# -----------------------------------------------------------------------------

@testset "Parrilo Example 3.1 --- " begin
	@polyvar x y
	p  = 2 * x * y^2 - x^2 - y
	Sx = @set 1 - x^2 >= 0
	Sy = @set 1 - y^2 >= 0

	expected_measure_x = AtomicMeasure(variables(Sx), [
		WeightedDiracMeasure([0.3968], 1)])
	expected_measure_y = AtomicMeasure(variables(Sy), [
		WeightedDiracMeasure([0.6299], 1)])
	expected_value = -0.4724

	value, measure_x, measure_y = solve_game(p, Sx, Sy, opt)

	eps = 1e-3
	@test isapprox(expected_value, value, atol=eps)
	@test isapprox(expected_measure_x, measure_x, atol=eps)
	@test isapprox(expected_measure_y, measure_y, atol=eps)
end

@testset "Parrilo Example 3.2 --- " begin
	@polyvar x y
	p  = 5x * y - 2x^2 - 2x * y^2 - y
	Sx = @set 1 - x^2 >= 0
	Sy = @set 1 - y^2 >= 0

	expected_measure_x = AtomicMeasure(variables(Sx), [
		WeightedDiracMeasure([0.2], 1)])
	expected_measure_y = AtomicMeasure(variables(Sy), [
		WeightedDiracMeasure([ 1], 0.78),
		WeightedDiracMeasure([-1], 0.22)])
	expected_value = -0.48

	value, measure_x, measure_y = solve_game(p, Sx, Sy, opt, iteration=0x1)

	eps = 1e-3
	@test isapprox(expected_value, value, atol=eps)
	@test isapprox(expected_measure_x, measure_x, atol=eps)
	@test isapprox(expected_measure_y, measure_y, atol=eps)
end
