using DataFrames
using CSV
using Plots
GR.inline("png")

function plot_curves(curve::Curve, tlabel::String, xlabel::String, ylabel::String, name::String)
    p = plot(curve.tvalues, curve.xvalues, label = "pressure", lc = :red, legend = :bottomright, xaxis = tlabel, yaxis = xlabel, show = false)
    plot!(twinx(), curve.tvalues, curve.yvalues, yaxis = ylabel, show = false, linestyle = :dash, label = false)
    plot!([-1], [0], xlims = extrema(curve.tvalues), lc = :blue, label = "mass") #evil hack from a random forum
    savefig(p, name)
end

using DataFrames
using CSV

function write_data(curve::Curve)
    df = DataFrame()
    df.radius = curve.tvalues
    df.pressure = curve.xvalues
    df.mass = curve.yvalues
    CSV.write("tov_data.csv", df)
end

function plot_from_datafile()
    df = CSV.File("tov_data.csv") |> DataFrame
    curve = Curve(df.radius, df.pressure, df.mass)
    plot_curves(curve)
end
