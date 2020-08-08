using DynamicPolynomials
using MultivariateMoments
using SemialgebraicSets
using CSDP

# -----------------------------------------------------------------------------
# Jiawang Nie, Zi Yang, & Guangming Zhou. (2018).
# The Saddle Point Problem of Polynomials.
# -----------------------------------------------------------------------------

# Players are inverted!

@testset "Saddle Example 6.1 (II) " begin
        @polyvar x[1:3] y[1:3]
        p =
                x[1]^3 + x[2]^3 - x[3]^3 - y[1]^3 - y[2]^3 +
                y[3]^3 +
                x[3]y[1]y[2] * (y[1] + y[2]) +
                x[2]y[1]y[3] * (y[1] + y[3]) +
                x[1]y[2]y[3] * (y[2] + y[3])
        Sx = @set sum(x) == 1 && x[1] >= 0 && x[2] >= 0 && x[3] >= 0
        Sy = @set sum(y) == 1 && y[1] >= 0 && y[2] >= 0 && y[3] >= 0
        opt = CSDP.Optimizer

        expected_measure_x = AtomicMeasure(
                variables(Sx),
                [WeightedDiracMeasure([0, 0, 1], 1)],
        )
        expected_measure_y = AtomicMeasure(
                variables(Sy),
                [WeightedDiracMeasure([0, 0, 1], 1)],
        )
        expected_value = 0

        pvalue, pmeasure_y, pmeasure_x =
                solve_game(+p, Sy, Sx, opt, iteration = 0x0)
        dvalue, dmeasure_x, dmeasure_y =
                solve_game(-p, Sx, Sy, opt, iteration = 0x0)


        eps = 1e-3
        @test isapprox(expected_value, pvalue, atol = eps)
        @test isapprox(dvalue, -pvalue, atol = eps)
        @test isapprox(expected_measure_y, pmeasure_y, atol = eps)
        @test isapprox(expected_measure_x, dmeasure_x, atol = eps)
end


@testset "Saddle Example 6.1 (III)" begin
        @polyvar x[1:4] y[1:4]
        p =
                sum([(x[i]^2 * y[i]^2) for i in 1:4]) - sum([
                        (i != j) ? (x[i]x[j] + y[i]y[j]) : 0
                        for i in 1:4, j in 1:4
                ])
        Sx = @set sum(x) == 1 &&
             x[1] >= 0 &&
             x[2] >= 0 &&
             x[3] >= 0 &&
             x[4] >= 0
        Sy = @set sum(y) == 1 &&
             y[1] >= 0 &&
             y[2] >= 0 &&
             y[3] >= 0 &&
             y[4] >= 0
        opt = CSDP.Optimizer

        expected_measure_x = AtomicMeasure(
                variables(Sx),
                [WeightedDiracMeasure([0.25, 0.25, 0.25, 0.25], 1)],
        )
        expected_measure_y = AtomicMeasure(
                variables(Sy),
                [
                        WeightedDiracMeasure([1, 0, 0, 0], 0.25),
                        WeightedDiracMeasure([0, 1, 0, 0], 0.25),
                        WeightedDiracMeasure([0, 0, 1, 0], 0.25),
                        WeightedDiracMeasure([0, 0, 0, 1], 0.25),
                ],
        )
        expected_value = -0.6875

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


@testset "Saddle Example 6.1 (IV) " begin
        @polyvar x[1:3] y[1:3]
        p =
                x[1]x[2]y[1]y[2] + x[2]x[3]y[2]y[3] + x[3]x[1]y[3]y[1] -
                x[1]^2 * y[3]^2 - x[2]^2 * y[1]^2 - x[3]^2 * y[2]^2
        Sx = @set sum(x) == 1 && x[1] >= 0 && x[2] >= 0 && x[3] >= 0
        Sy = @set sum(y) == 1 && y[1] >= 0 && y[2] >= 0 && y[3] >= 0
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
                [WeightedDiracMeasure([1 / 3, 1 / 3, 1 / 3], 1)],
        )
        expected_value = -1 / 9

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

@testset "Saddle Example 6.3 (I)  " begin
        @polyvar x[1:3] y[1:3]
        p =
                sum([(x[i] + y[i]) for i in 1:3]) - prod([(x[i] - y[i]) for i in 1:3])
        Sx = @set x[1]^2 <= 1 && x[2]^2 <= 1 && x[3]^2 <= 1
        Sy = @set y[1]^2 <= 1 && y[2]^2 <= 1 && y[3]^2 <= 1
        opt = CSDP.Optimizer

        expected_measure_x = AtomicMeasure(
                variables(Sx),
                [
                        WeightedDiracMeasure([-1, -1, 1], 1 / 3),
                        WeightedDiracMeasure([1, -1, -1], 1 / 3),
                        WeightedDiracMeasure([-1, 1, -1], 1 / 3),
                ],
        )
        expected_measure_y = AtomicMeasure(
                variables(Sy),
                [WeightedDiracMeasure([1, 1, 1], 1)],
        )
        expected_value = 2

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

@testset "Saddle Example 6.3 (II)  " begin
        @polyvar x[1:3] y[1:3]
        p =
                y'y - x'x +
                sum([(x[i]y[j] - x[j]y[i]) for i = 1:3, j = 1:3 if i < j])
        Sx = @set x[1]^2 <= 1 && x[2]^2 <= 1 && x[3]^2 <= 1
        Sy = @set y[1]^2 <= 1 && y[2]^2 <= 1 && y[3]^2 <= 1
        opt = CSDP.Optimizer

        expected_measure_x = AtomicMeasure(
                variables(Sx),
                [
                        WeightedDiracMeasure([-1.0, -1.0, -1.0], 1 / 8),
                        WeightedDiracMeasure([1.0, 1.0, 1.0], 1 / 8),
                        WeightedDiracMeasure([1.0, -1.0, -1.0], 1 / 8),
                        WeightedDiracMeasure([-1.0, 1.0, -1.0], 1 / 8),
                        WeightedDiracMeasure([-1.0, -1.0, 1.0], 1 / 8),
                        WeightedDiracMeasure([-1.0, 1.0, 1.0], 1 / 8),
                        WeightedDiracMeasure([1.0, 1.0, -1.0], 1 / 8),
                        WeightedDiracMeasure([1.0, -1.0, 1.0], 1 / 8),
                ],
        )
        expected_measure_y = AtomicMeasure(
                variables(Sy),
                [
                        WeightedDiracMeasure([1.0, 1.0, 1.0], 1 / 8),
                        WeightedDiracMeasure([1.0, 1.0, -1.0], 1 / 8),
                        WeightedDiracMeasure([-1.0, 1.0, 1.0], 1 / 8),
                        WeightedDiracMeasure([-1.0, 1.0, -1.0], 1 / 8),
                        WeightedDiracMeasure([-1.0, -1.0, -1.0], 1 / 8),
                        WeightedDiracMeasure([1.0, -1.0, 1.0], 1 / 8),
                        WeightedDiracMeasure([-1.0, -1.0, 1.0], 1 / 8),
                        WeightedDiracMeasure([1.0, -1.0, -1.0], 1 / 8),
                ],
        )
        expected_value = 0

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



@testset "Saddle Example 6.4 (I)  " begin
        @polyvar x[1:3] y[1:3]
        p =
                x[1]^3 +
                x[2]^3 +
                x[3]^3 +
                y[1]^3 +
                y[2]^3 +
                y[3]^3 +
                2 * (x[1]x[2]y[1]y[2] + x[1]x[3]y[1]y[3] + x[2]x[3]y[2]y[3])
        Sx = @set x'x == 1
        Sy = @set y'y == 1
        opt = CSDP.Optimizer

        expected_measure_x = AtomicMeasure(
                variables(Sx),
                [
                        WeightedDiracMeasure([0, -1, 0], 1 / 3),
                        WeightedDiracMeasure([-1, 0, 0], 1 / 3),
                        WeightedDiracMeasure([0, 0, -1], 1 / 3),
                ],
        )
        expected_measure_y = AtomicMeasure(
                variables(Sy),
                [
                        WeightedDiracMeasure([0, 1, 0], 1 / 3),
                        WeightedDiracMeasure([1, 0, 0], 1 / 3),
                        WeightedDiracMeasure([0, 0, 1], 1 / 3),
                ],
        )
        expected_value = 0

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

@testset "Saddle Example 6.5 (I)  " begin
        @polyvar x[1:3] y[1:3]
        p =
                x[1]^2 * y[1] + 2 * x[2]^2 * y[2] + 3 * x[3]^2 * y[3] - x[1] -
                x[2] - x[3]
        Sx = @set x'x <= 1
        Sy = @set y'y <= 1
        opt = CSDP.Optimizer

        expected_measure_x = AtomicMeasure(
                variables(Sx),
                [WeightedDiracMeasure([0.726, 0.458, 0.349], 1)],
        )
        expected_measure_y = AtomicMeasure(
                variables(Sy),
                [WeightedDiracMeasure([0.688, 0.546, 0.477], 1)],
        )
        expected_value = -0.7666

        pvalue, pmeasure_y, pmeasure_x =
                solve_game(+p, Sy, Sx, opt, iteration = 0x0)
        dvalue, dmeasure_x, dmeasure_y =
                solve_game(-p, Sx, Sy, opt, iteration = 0x0)


        eps = 1e-3
        @test isapprox(expected_value, pvalue, atol = eps)
        @test isapprox(dvalue, -pvalue, atol = eps)
        @test isapprox(expected_measure_y, pmeasure_y, atol = eps)
        @test isapprox(expected_measure_x, dmeasure_x, atol = eps)
end
