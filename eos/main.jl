include("../constants.jl")
include("./polytropic.jl")

using Plots
function plot_pϵ(pmin::Real, pmax::Real)
    p = range(pmin, pmax, length = 1000)
    ϵ = broadcast(polytrope, p, γ_nonrel, K_NONREL)

    pl = plot(ϵ, p)
    xlabel!(pl, raw"$\epsilon$", dpi = 600)
    ylabel!(pl, raw"$p$")
    savefig(pl, "eos_plot.png")
end
