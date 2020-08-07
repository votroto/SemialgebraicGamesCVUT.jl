# SemialgebraicGames.jl
Solver for two-player zero-sum polynomial games on semialgebraic sets.
Contrary to the repository name, this package is called `SemialgebraicGames`.

> Tip: Julia is "fast" is its own _unique_ way, though you may want to keep a book nearby, just in case...
## To add
Open your _Julia_ interpreter, press `]` to enter the _pkg_-mode, and run
```
add https://github.com/votroto/SemialgebraicGamesCVUT.jl.git
```
to add this package (or use `dev` instead of `add` to keep the package "alive" for development).

Finaly, press `backspace` to exit the _pkg_-mode.

## To use
```julia
# Load this package. (This will take a while)
using SemialgebraicGames
# Packages to define the payoff function, and the strategy sets.
using DynamicPolynomials
using SemialgebraicSets
# Example semidefinite optimizer.
using CSDP
opt = CSDP.Optimizer

# Example game definition
@polyvar x y
p  = (x - y)^2
Sx = @set 1 - x^2 >= 0
Sy = @set 1 - y^2 >= 0

# Run the solver. (First time load. Again.)
value_x, measure_x, measure_y = solve_game(p, Sx, Sy, opt, iteration=0x1)
```

## To remove:
Open your _Julia_ interpreter, press `]` to enter the _pkg_-mode, and run
```
rm SemialgebraicGames
```
Finaly, press `backspace` to exit the _pkg_-mode.