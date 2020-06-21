using DynamicPolynomials
using MultivariateMoments
using SemialgebraicSets
using SCS

# -----------------------------------------------------------------------------
# Trivial tests
# -----------------------------------------------------------------------------

@testset "Guessing game --------- " begin
	@polyvar x y
	p  = (x - y)^2
	Sx = @set 1 - x^2 >= 0
	Sy = @set 1 - y^2 >= 0

	expected_measure_x = AtomicMeasure(variables(Sx), [
		WeightedDiracMeasure([1], 0.5),
		WeightedDiracMeasure([-1], 0.5)])
	expected_measure_y = AtomicMeasure(variables(Sy), [
		WeightedDiracMeasure([0], 1)])
	expected_value_y = 1

	actual_value_x, actual_measure_x = solve_game(-p, Sy, Sx, SCS.Optimizer, iteration=1)
	actual_value_y, actual_measure_y = solve_game(p, Sx, Sy, SCS.Optimizer, iteration=1)

	eps = 1e-3
	@test isapprox(actual_value_x, -actual_value_y, atol=eps)
	@test isapprox(expected_value_y, actual_value_y, atol=eps)
	@test isapprox(expected_measure_y, actual_measure_y, atol=eps)
	@test isapprox(expected_measure_x, actual_measure_x, atol=eps)
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
	expected_value_y = -0.4724

	actual_value_x, actual_measure_x = solve_game(-p, Sy, Sx, SCS.Optimizer)
	actual_value_y, actual_measure_y = solve_game(p, Sx, Sy, SCS.Optimizer)

	eps = 1e-3
	@test isapprox(actual_value_x, -actual_value_y, atol=eps)
	@test isapprox(expected_value_y, actual_value_y, atol=eps)
	@test isapprox(expected_measure_y, actual_measure_y, atol=eps)
	@test isapprox(expected_measure_x, actual_measure_x, atol=eps)
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
	expected_value_y = -0.48

	actual_value_x, actual_measure_x = solve_game(-p, Sy, Sx, SCS.Optimizer, iteration=1)
	actual_value_y, actual_measure_y = solve_game(p, Sx, Sy, SCS.Optimizer, iteration=1)

	eps = 1e-3
	@test isapprox(actual_value_x, -actual_value_y, atol=eps)
	@test isapprox(expected_value_y, actual_value_y, atol=eps)
	@test isapprox(expected_measure_y, actual_measure_y, atol=eps)
	@test isapprox(expected_measure_x, actual_measure_x, atol=eps)
end

# -----------------------------------------------------------------------------
# Jiawang Nie, Zi Yang, & Guangming Zhou. (2018).
# The Saddle Point Problem of Polynomials.
# -----------------------------------------------------------------------------

@testset "Saddle Example 6.1 (I)  " begin
	@polyvar x[1:3] y[1:3]
	p = x[1]x[2] + x[2]x[3] + x[3]y[1] + x[1]y[3] + y[1]y[2] + y[2]y[3]
	Sx = @set sum(x) == 1 && x[1] >= 0 && x[2] >= 0 && x[3] >= 0
	Sy = @set sum(y) == 1 && y[1] >= 0 && y[2] >= 0 && y[3] >= 0

	expected_measure_x = AtomicMeasure(variables(Sy), [
		WeightedDiracMeasure([0.25,0.5,0.25], 1)])
	expected_measure_y = AtomicMeasure(variables(Sx), [
		WeightedDiracMeasure([0,1,0], 1)])
	expected_value_y = 0.25

	actual_value_x, actual_measure_x = solve_game(-p, Sy, Sx, SCS.Optimizer, iteration=1)
	actual_value_y, actual_measure_y = solve_game(p, Sx, Sy, SCS.Optimizer, iteration=1)

	eps = 1e-3
	@test isapprox(actual_value_x, -actual_value_y, atol=eps)
	@test isapprox(expected_value_y, actual_value_y, atol=eps)
	#@test isapprox(expected_measure_y, actual_measure_y, atol = eps)
	#@test isapprox(expected_measure_x, actual_measure_x, atol = eps)
end

@testset "Saddle Example 6.1 (II) " begin
	@polyvar x[1:3] y[1:3]
	p = x[1]^3 + x[2]^3 - x[3]^3 - y[1]^3 - y[2]^3 + y[3]^3
		+ x[3]y[1]y[2]*(y[1] + y[2]) + x[2]y[1]y[3]*(y[1] + y[3])
		+ x[1]y[2]y[3]*(y[2] + y[3])
	Sx = @set sum(x) == 1 && x[1] >= 0 && x[2] >= 0 && x[3] >= 0
	Sy = @set sum(y) == 1 && y[1] >= 0 && y[2] >= 0 && y[3] >= 0

	expected_measure_x = AtomicMeasure(variables(Sx), [
		WeightedDiracMeasure([0,0,1], 1)])
	expected_measure_y = AtomicMeasure(variables(Sy), [
		WeightedDiracMeasure([0,0,1], 1)])
	expected_value_y = 0

	actual_value_x, actual_measure_x = solve_game(-p, Sy, Sx, SCS.Optimizer, iteration=1)
	actual_value_y, actual_measure_y = solve_game(p, Sx, Sy, SCS.Optimizer, iteration=1)

	eps = 1e-3
	@test isapprox(actual_value_x, -actual_value_y, atol=eps)
	@test isapprox(expected_value_y, actual_value_y, atol=eps)
	#@test isapprox(expected_measure_y, actual_measure_y, atol = eps)
	#@test isapprox(expected_measure_x, actual_measure_x, atol = eps)
end

@testset "Saddle Example 6.1 (III)" begin
	@polyvar x[1:4] y[1:4]
	p = sum([(x[i]^2*y[i]^2) for i in 1:4])	- sum([(i!=j) ? (x[i]x[j] + y[i]y[j]) : 0 for i=1:4, j=1:4])
	Sx = @set sum(x) == 1 && x[1] >= 0 && x[2] >= 0 && x[3] >= 0 && x[4] >= 0
	Sy = @set sum(y) == 1 && y[1] >= 0 && y[2] >= 0 && y[3] >= 0 && y[4] >= 0

	expected_measure_x = AtomicMeasure(variables(Sx), [
		WeightedDiracMeasure([0.25,0.25,0.25,0.25], 1)])
	expected_measure_y = AtomicMeasure(variables(Sy), [
		WeightedDiracMeasure([1,0,0,0], 0.25),
		WeightedDiracMeasure([0,1,0,0], 0.25),
		WeightedDiracMeasure([0,0,1,0], 0.25),
		WeightedDiracMeasure([0,0,0,1], 0.25)])
	expected_value_y = -0.6875

	actual_value_x, actual_measure_x = solve_game(-p, Sy, Sx, SCS.Optimizer, iteration=1)
	actual_value_y, actual_measure_y = solve_game(p, Sx, Sy, SCS.Optimizer, iteration=1)

	eps = 1e-3
	@test isapprox(actual_value_x, -actual_value_y, atol=eps)
	@test isapprox(expected_value_y, actual_value_y, atol=eps)
	#@test isapprox(expected_measure_y, actual_measure_y, atol = eps)
	#@test isapprox(expected_measure_x, actual_measure_x, atol = eps)
end

@testset "Saddle Example 6.2 (I)  " begin
	@polyvar x[1:2] y[1:2]
	p = (x[1] + x[2] + y[1] + y[2] + 1)^2 - 4 * (x[1]x[2] + x[2]y[1] + y[1]y[2] + y[2] + x[1])
	Sx = @set -x[1]^2 + x[1] >= 0 && -x[2]^2 + x[2] >= 0
	Sy = @set -y[1]^2 + y[1] >= 0 && -y[2]^2 + y[2] >= 0

	expected_measure_x = AtomicMeasure(variables(Sx), [
		WeightedDiracMeasure([0.3249, 0.3249], 1)])
	expected_measure_y = AtomicMeasure(variables(Sy), [
		WeightedDiracMeasure([1, 0], 1)])
	expected_value_y = 4

	actual_value_x, actual_measure_x = solve_game(-p, Sy, Sx, SCS.Optimizer, iteration=1)
	actual_value_y, actual_measure_y = solve_game(p, Sx, Sy, SCS.Optimizer, iteration=1)

	eps = 1e-3
	@test isapprox(actual_value_x, -actual_value_y, atol=eps)
	@test isapprox(expected_value_y, actual_value_y, atol=eps)
	#@test isapprox(expected_measure_y, actual_measure_y, atol = eps)
	#@test isapprox(expected_measure_x, actual_measure_x, atol = eps)
end

@testset "Saddle Example 6.3 (I)  " begin
	@polyvar x[1:3] y[1:3]
	p = sum([(x[i] + y[i]) for i=1:3]) - prod([(x[i] - y[i]) for i=1:3])
	Sx = @set x[1]^2 <= 1 && x[2]^2 <= 1 && x[3]^2 <= 1
	Sy = @set y[1]^2 <= 1 && y[2]^2 <= 1 && y[3]^2 <= 1

	expected_measure_x = AtomicMeasure(variables(Sx), [
		WeightedDiracMeasure([-1,-1, 1], 1/3),
		WeightedDiracMeasure([-1, 1, 1], 1/3),
		WeightedDiracMeasure([ 1,-1,-1], 1/3)])
	expected_measure_y = AtomicMeasure(variables(Sy), [
		WeightedDiracMeasure([ 1, 1, 1], 1)])
	expected_value_y = -2

	actual_value_x, actual_measure_x = solve_game(-p, Sy, Sx, SCS.Optimizer, iteration=1)
	actual_value_y, actual_measure_y = solve_game(p, Sx, Sy, SCS.Optimizer, iteration=1)

	eps = 1e-3
	@test isapprox(actual_value_x, -actual_value_y, atol=eps)
	@test isapprox(expected_value_y, actual_value_y, atol=eps)
end

@testset "Saddle Example 6.5 (I)  " begin
	@polyvar x[1:3] y[1:3]
	p = x[1]^2*y[1] + 2*x[2]^2*y[2] + 3*x[3]^2*y[3] - x[1] - x[2] - x[3]
	Sx = @set x'x <= 1
	Sy = @set y'y <= 1

	expected_value_y = 0.7666

	actual_value_x, actual_measure_x = solve_game(-p, Sy, Sx, SCS.Optimizer, iteration=1)
	actual_value_y, actual_measure_y = solve_game(p, Sx, Sy, SCS.Optimizer, iteration=1)

	eps = 1e-3
	@test isapprox(actual_value_x, -actual_value_y, atol=eps)
	@test isapprox(expected_value_y, actual_value_y, atol=eps)
end
