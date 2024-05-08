@testset "Differential Equation Solver" begin

α = 0.01
β = 0.02

f(t, x) = [(1 - α*x[2])*x[1], (β*x[1] - 1)*x[2]]
sol = solvesystem(f, 0.0, [20.0, 20.0], 0.01, (i, t, x) -> t <= 15)
plot(sol[:, 1], sol[:, 2], label="x1", xlabel="t")
plot!(sol[:, 1], sol[:, 3], label="x2", xlabel="t")
savefig(joinpath(dirname(@__FILE__), "out", "test_diffsolve.png"))

end