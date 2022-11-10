include("diff.jl") #useless include 
include("tov.jl")

using Plots
using DataFrames
using CSV

function plot_functions(curve::Curve)
    rticks = [5000, 10000, 15000, 20000, 25000, 30000]

    p = plot(curve.tvalues, curve.xvalues, label = false, xticks = rticks)
    xlabel!(p, raw"$r (km)$")
    ylabel!(p, raw"$p (J/m^3)$")
    savefig(p, "pressure_plot.png")

    m = plot(curve.tvalues, curve.yvalues, label = false, xticks = rticks)
    xlabel!(m, raw"$r (km)$")
    ylabel!(m, raw"$mass (M⊙)$")
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

function solve(p₀::Real)::Curve
    curve = try solve_tov(p₀)
        catch err
            println(err)
            return
        end

    write_data(curve)

    return curve
end

function solve_plot(p₀::Real)
    curve = solve(p₀)

    plot_functions(curve)
end

function solve_data()
   solve()
   return
end

using Interpolations

function solve_star_curve(pa::Real, pb::Real)
    n = 10000
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
