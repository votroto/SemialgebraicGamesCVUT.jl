using DynamicPolynomials
using MultivariateMoments
using SemialgebraicSets
using CSDP

# -----------------------------------------------------------------------------
# Jiawang Nie, Zi Yang, & Guangming Zhou. (2018).
# The Saddle Point Problem of Polynomials.
# -----------------------------------------------------------------------------

@testset "Bad Saddle Example 6.1 (I)" begin
        @polyvar x[1:3] y[1:3]
        p = x[1]x[2] + x[2]x[3] + x[3]y[1] + x[1]y[3] + y[1]y[2] + y[2]y[3]
        Sx = @set sum(x) == 1 && x[1] >= 0 && x[2] >= 0 && x[3] >= 0
        Sy = @set sum(y) == 1 && y[1] >= 0 && y[2] >= 0 && y[3] >= 0
        opt = CSDP.Optimizer

        expected_measure_x = AtomicMeasure(
                variables(Sx),
                [WeightedDiracMeasure([0, 1, 0], 1)],
        )
        expected_measure_y = AtomicMeasure(
                variables(Sy),
                [
                        WeightedDiracMeasure([0.02, 0.5, 0.48], 0.3),
                        WeightedDiracMeasure([0.25, 0.5, 0.25], 0.3),
                        WeightedDiracMeasure([0.48, 0.5, 0.02], 0.3),
                ],
        )
        expected_value = 0.25

        pvalue, pmeasure_y, pmeasure_x =
                solve_game(+p, Sy, Sx, opt, iteration = 0x2)
        dvalue, dmeasure_x, dmeasure_y =
                solve_game(-p, Sx, Sy, opt, iteration = 0x1)


        eps = 1e-3
        @test isapprox(expected_value, pvalue, atol = eps)
        @test isapprox(dvalue, -pvalue, atol = eps)
        @test isapprox(expected_measure_y, pmeasure_y, atol = 1e-1) # WARNING!
        @test isapprox(expected_measure_x, dmeasure_x, atol = eps)
end


@testset "Bad Saddle Example 6.2 (I)  " begin
        @polyvar x[1:2] y[1:2]
        p =
                (x[1] + x[2] + y[1] + y[2] + 1)^2 -
                4 * (x[1]x[2] + x[2]y[1] + y[1]y[2] + y[2] + x[1])
        Sx = @set x[1] >= 0 && x[2] >= 0 && 1 - x[1] >= 0 && 1 - x[2] >= 0
        Sy = @set y[1] >= 0 && y[2] >= 0 && 1 - y[1] >= 0 && 1 - y[2] >= 0
        opt = CSDP.Optimizer

        expected_measure_y =
                AtomicMeasure(variables(Sy), [WeightedDiracMeasure([1, 0], 1)])
        expected_value = 4

        pvalue, pmeasure_y, pmeasure_x =
                solve_game(+p, Sy, Sx, opt, iteration = 0x1)
        dvalue, dmeasure_x, dmeasure_y =
                solve_game(-p, Sx, Sy, opt, iteration = 0x1)


        eps = 1e-3
        @test isapprox(expected_value, pvalue, atol = eps)
        @test isapprox(dvalue, -pvalue, atol = eps)
        @test isapprox(expected_measure_y, pmeasure_y, atol = eps)
end

@testset "Unknown Saddle Example 6.2 (II)  " begin
        @polyvar x[1:3] y[1:3]
        p =
                sum([(x[i] + y[i]) for i in 1:3]) - prod([(x[i] - y[i]) for i in 1:3])
        Sx = @set x[1] >= 0 &&
             x[2] >= 0 &&
             x[3] >= 0 &&
             1 - x[1] >= 0 &&
             1 - x[2] >= 0 &&
             1 - x[3] >= 0
        Sy = @set y[1] >= 0 &&
             y[2] >= 0 &&
             y[3] >= 0 &&
             1 - y[1] >= 0 &&
             1 - y[2] >= 0 &&
             1 - y[3] >= 0
        opt = CSDP.Optimizer

        expected_measure_x = AtomicMeasure(
                variables(Sx),
                [
                        WeightedDiracMeasure([0, 1, 0], 1 / 3),
                        WeightedDiracMeasure([1, 0, 0], 1 / 3),
                        WeightedDiracMeasure([0, 0, 1], 1 / 3),
                ],
        )
        expected_measure_y = AtomicMeasure(
                variables(Sy),
                [WeightedDiracMeasure([1, 1, 1], 1)],
        )
        expected_value = 4

        pvalue, pmeasure_y, pmeasure_x =
                solve_game(+p, Sy, Sx, opt, iteration = 0x1)
        dvalue, dmeasure_x, dmeasure_y =
                solve_game(-p, Sx, Sy, opt, iteration = 0x1)


        eps = 1e-3
        @test isapprox(expected_value, pvalue, atol = eps)
        @test isapprox(dvalue, -pvalue, atol = eps)
        @test isapprox(expected_measure_y, pmeasure_y, atol = eps)
        @test isapprox(expected_measure_x, dmeasure_x, atol = eps)
end


#=
#@testset "Saddle Example 6.4 (II)  " begin
#	no saddle point - useless test
#	@polyvar x[1:3] y[1:3]
#	p = x[1]^2y[1]^2 + x[2]^2y[2]^2 + x[3]^2y[3]^2 + x[1]^2y[2]y[3] + x[2]^2y[1]y[3] + x[3]^2y[1]y[2] + y[1]^2x[2]x[3] + y[2]^2x[1]x[3] + y[3]^2x[1]x[2]
#	Sx = @set x'x == 1
#	Sy = @set y'y == 1
#       opt = CSDP.Optimizer
#
#	expected_value_y = 0.7666
#
#	actual_value_x, actual_measure_x = solve_game(-p, Sy, Sx, opt, iteration=1)
#	actual_value_y, actual_measure_y = solve_game(p, Sx, Sy, opt, iteration=1)
#
#	eps = 1e-3
#	@test isapprox(actual_value_x, -actual_value_y, atol=eps)
#	@test isapprox(expected_value_y, actual_value_y, atol=eps)
#end

=#
