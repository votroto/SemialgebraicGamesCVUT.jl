function gale_gross(x, y, spectrum_x, spectrum_y, weights_x, weights_y)
	function _f(i, x, sx)
		if i == 0
			prod([(x - s)^2 for s in sx])
		else
			prod([(x - s)^2 / (sx[i] - s)^2 for s in sx if s != sx[i]])
		end
	end

	function _ϕ(x, as)
		prod([(x - s)^2 for s in as])
	end

	function _ϕ(i, x, as)
		prod([(x - s)^2 for s in as if s != as[i + 1]])
	end

	μs = weights_x
	νs = weights_y
	m = length(μs)
	n = length(νs)

	αs = rand(Float64, n + 1) * 2 .- 1
	βs = rand(Float64, m + 1) * 2 .- 1

	f(i, x) = _f(i, x, spectrum_x)
	ϕ(i, x) = _ϕ(i, x, αs)
	ϕ(x) = _ϕ(x, αs)
	g(i, y) = _f(i, y, spectrum_y)
	ψ(i, y) = _ϕ(i, y, βs)
	ψ(y) = _ϕ(y, βs)


	m1 = +f(0, x)ϕ(x) * (g(0, y)ϕ(0, x) + sum([(g(j, y) - νs[j]) * ϕ(j, x) for j = 1:n]))
	m2 = -g(0, y)ψ(y) * (f(0, x)ψ(0, y) + sum([(f(i, x) - μs[i]) * ψ(i, y) for i = 1:m]))
	m3 = -(f(0, x)ϕ(x))^2 + (g(0, y)ψ(y))^2

	m1 + m2 + m3
end