include("./polytropic.jl")

#TODO: make dimensionless  ̄p and  ̄ϵ and choose appropriate values of α and β

using Plots
function plot_pϵ()
    p = range(1.0e13, 1.0e-22, length = 100)
    ϵ = broadcast(ϵ_rel, p)

    pl = plot(ϵ, p)
    xlabel!(pl, raw"$\epsilon$", dpi = 600)
    ylabel!(pl, raw"$p$")
    savefig(pl, "eos_plot.png")
end
