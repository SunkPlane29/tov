include("diff.jl") #useless include
include("tov.jl")
include("util.jl")

using Plots
using DataFrames
using CSV

function plot_functions(pcurve::Curve, mcurve::Curve)
    p = plot(pcurve.tvalues, pcurve.yvalues, label="pressure")
    xlabel!(p, "r in km")
    ylabel!(p, "p in fm^(-4)")
    savefig(p, "pressure_plot.png")

    m = plot(mcurve.tvalues, mcurve.yvalues, label="mass")
    xlabel!(m, "r in km")
    ylabel!(m, "mass in M⊙")
    savefig(m, "mass_plot.png")
end

function write_data(pcurve::Curve, mcurve::Curve)
    df = DataFrame()
    df.radius = pcurve.tvalues
    df.pressure = pcurve.yvalues
    df.mass = mcurve.yvalues
    CSV.write("tov_data.csv", df)
end

function plot_from_datafile()
    df = CSV.File("tov_data.csv") |> DataFrame
    pcurve = Curve(df.radius, df.pressure)
    mcurve = Curve(df.radius, df.mass)
    plot_functions(pcurve, mcurve)
end

function solve_write()::Tuple{Curve, Curve}
    (pcurve, mcurve) = try solve_tov()
        catch err
            println(err)
            return
            end

    write_data(pcurve, mcurve)

    return (pcurve, mcurve)
end

function solve_normal()
    (pcurve, mcurve) = solve_write()

    plot_functions(pcurve, mcurve)
end

function solve_data()
   solve_write()
   return
end

function main()
    if get_only_graph()
        plot_from_datafile()
        return
    end

    if get_only_data() solve_data(); return end
    solve_normal()
end

main()
