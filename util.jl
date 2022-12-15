using Plots
GR.inline("png")

function plot_curves(curve::Curve, tlabel::String, xlabel::String, ylabel::String, xname::String, yname::String)
    p = plot(curve.tvalues, curve.xvalues, label = false, show = false)
    xlabel!(p, tlabel)
    ylabel!(p, xlabel)
    savefig(p, xname)

    m = plot(curve.tvalues, curve.yvalues, label = false, show = false)
    xlabel!(m, tlabel)
    ylabel!(m, ylabel)
    savefig(m, yname)
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

macro run(pressure, relativistic = true)
    return :(include("main.jl") ; solve_plot($pressure, $relativistic))
end
