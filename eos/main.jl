include("../constants.jl")
include("./polytropic.jl")

using Plots

function plot_pϵ(pmin::Real, pmax::Real)
    p = range(pmin, pmax, length = 1000)
    ϵ = broadcast(polytrope, p, γ_nonrel, K_NONREL)

    pl = plot(ϵ, p, label = "p", legendposition = :bottomright)
    plot!(pl, p, ϵ, label = "ϵ")

    pt = (0.752, 0.752)
    plot!(pl, [pt[1]], [pt[2]], markershape = :circle, markercolor = :black, markerstrokewidth = 0, markersize = 2)

    ptstr = string(pt)
    xvalstr = string(pt[1])
    annotate!(pl, pt..., text(xvalstr, 14, :bottom, :left, :black))

    savefig(pl, "eos_plot.png")
end
