include("diff.jl") #useless include
include("tov.jl")

using Plots

function main()
    (pcurve, mcurve) = try solve_tov()
        catch err
            println(err)
            return
        end

    p = plot(pcurve.tvalues, pcurve.yvalues, label="pressure")
    xlabel!(p, "r in km")
    ylabel!(p, "p in fm^(-4)")
    savefig(p, "plot.png")
end

main()
