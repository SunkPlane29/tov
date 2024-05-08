@testset "Cubic Spline Interpolation" begin

x = [-1.0, 0.0, 3.0]
y = [0.5, 0.0, 3.0]

cs = CubicSplineInterpolation(x, y)

plot(x, y, seriestype=:scatter, label="data", xlabel="x", ylabel="y")

xinterp = range(x[1], x[end], length=100)
plot!(xinterp, cs.(xinterp), label="cubic spline interpolation")
savefig(joinpath(dirname(@__FILE__), "out", "test_cubicspline.png"))

end