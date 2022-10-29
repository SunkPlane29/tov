include("diff.jl") #useless include
include("tov.jl")
include("util.jl")

using Plots
using DataFrames
using CSV

function plot_functions(curve::Curve)
    p = plot(curve.tvalues, curve.xvalues, label = false)
    xlabel!(p, "r in km")
    ylabel!(p, "p (dimensionless)")
    savefig(p, "pressure_plot.png")

    m = plot(curve.tvalues, curve.yvalues, label = false)
    xlabel!(m, "r in km")
    ylabel!(m, "mass in M⊙")
    savefig(m, "mass_plot.png")
end

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
    plot_functions(curve)
end

function solve_write()::Curve
    curve = try solve_tov(1.0e-15)
        catch err
            println(err)
            return
        end

    write_data(curve)

    return curve
end

function solve_normal()
    curve = solve_write()

    plot_functions(curve)
end

function solve_data()
   solve_write()
   return
end

function solve_star_curve()
    pa = 1e-15
    pb = 1e-13
    n = 100
    h = (pb - pa)/n

    Rvalues = []
    Mvalues = []

    for i = 1:n
        curve = try solve_tov(pa + (i-1)*h)
            catch err
                println(err)
                return
            end

        append!(Rvalues, last(curve.tvalues))
        append!(Mvalues, last(curve.yvalues))
    end

    p = plot(Rvalues, Mvalues, legend = false)
    xlabel!(p, "Radius (km)")
    ylabel!(p, "Mass (M⊙)")
    savefig("tov_plot.png")
end
