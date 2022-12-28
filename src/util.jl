using DataFrames
using CSV
using Plots
GR.inline("png")

#plot_curves plot the curve object (defined in diff.jl) with the required xlabel, ylabel, tlabel and name of the plot
#since both x and y curves depend on the same value of the tcurve, the curves are plotted using the same xaxis, but
#different yaxis
function plot_curves(curve::Curve, tlabel::String, xlabel::String, ylabel::String, name::String)
    p = plot(curve.tvalues, curve.xvalues, label = "pressure", lc = :red, legend = :bottomright, xaxis = tlabel, yaxis = xlabel, show = false)
    plot!(twinx(), curve.tvalues, curve.yvalues, yaxis = ylabel, show = false, linestyle = :dash, label = false)
    plot!([-1], [0], xlims = extrema(curve.tvalues), lc = :blue, label = "mass") #evil hack from a random forum
    savefig(p, name)
end

using DataFrames
using CSV

#TODO: there should be somewhere denoting the units
#write_data writes the curve data to a datafile, note that this function is not really general, so this will only work
#in the context of TOV
function write_data(curve::Curve)
    df = DataFrame()
    df.radius = curve.tvalues
    df.pressure = curve.xvalues
    df.mass = curve.yvalues
    CSV.write("single_star_data.csv", df)
end

#plot_from_datafile plots the TOV (single star curve) from a csv datafile
function plot_from_datafile()
    df = CSV.File("tov_data.csv") |> DataFrame
    curve = Curve(df.radius, df.pressure, df.mass)
    plot_curves(curve)
end