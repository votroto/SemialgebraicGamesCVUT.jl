# SemialgebraicGames.jl
`SemialgebraicGames.jl` is a solver for two-player zero-sum polynomial games on semialgebraic sets. It also includes DSOS and SDSOS optimization capabilities, along with an iterative change-of-basis algorithm. Contrary to the repository name, this package is called `SemialgebraicGames.jl`.

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
The problem can also be solved using an iterative change-of-basis technique based on DDP and SDDP.
```julia
using ECOS
opt = () -> CoBModel(ECOS.Optimizer, SDD, iterations=0x5)
```
> Keep in mind that SDDP and DDP are conservative relaxations. Therefore, `solve_game(p, Sx, Sy, ...)` is unlikely to equal the dual `solve_game(-p, Sy, Sx, ...)`. Try both.

It is also possible to use LP and SOCP solvers simply by imbueing them with diagonal-dominance bridging capabilities. This results in the standard DDP and SDDP optimization.
```julia
opt = lazy_relax(ECOS.Optimizer, SDDPBridge)
```

## To remove:
Open your _Julia_ interpreter, press `]` to enter the _pkg_-mode, and run
```
rm SemialgebraicGames
```
Finaly, press `backspace` to exit the _pkg_-mode.

## Acknowledgement
- Ahmadi, A., & Majumdar, A. (2017). Dsos and sdsos optimization: More tractable alternatives to sum of squares and semidefinite optimization. SIAM Journal on Applied Algebra and Geometry, 3. https://doi.org/10.1137/18M118935X
- Lasserre, J.-B. (2004). Global optimization with polynomials and the problem of moments. SIAM Journal on Optimization, 11. https://doi.org/10.1137/S1052623400366802
- Parrilo, P. (2007). Polynomial games and sum of squares optimization, 2855–2860. https://doi.org/10.1109/CDC.2006.377261

