include("tov.jl")
include("eos/polytropic.jl")

using Plots
GR.inline("png")
using DataFrames
using CSV

function plot_functions(curve::Curve)
    rticks = [5000, 10000, 15000, 20000, 25000, 30000]

    p = plot(curve.tvalues, curve.xvalues, label = false, xticks = rticks, show = false)
    xlabel!(p, raw"$r$ (km)")
    ylabel!(p, raw"$p$ (J/m$^3$)")
    savefig(p, "pressure_plot.png")

    m = plot(curve.tvalues, curve.yvalues, label = false, xticks = rticks, show = false)
    xlabel!(m, raw"$r$ (km)")
    ylabel!(m, raw"$M$ (M$_\odot$)")
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

function solve(p₀::Real, ϵ₀::Real, r₀::Real, eos::Function)::Curve
    curve = try solve_tov(p₀, ϵ₀, r₀, eos)
        catch err
            println(err)
            return
        end

    write_data(curve)

    return curve
end

#SUGESTION: use p₀ = 1.54e-16, ϵ₀ = 5.61970127e+38 and r₀ = 0.9319e4
#TODO: make some code later to pick these parameters to fit observations (? i don't know if it is allowed ?)
function solve_plot(p₀::Real, ϵ₀::Real, r₀::Real)
    #this make simpler to change from relativistic to non-relativistic later
    γ = γ_rel
    polytrope(p) = rel_polytrope(p)

    eos_const = 1/ϵ₀^((γ-1)/γ)
    eos(p) = p <= 0 ? 0 : polytrope(p)*eos_const
    curve = solve(p₀, ϵ₀, r₀, eos)

    plot_functions(curve)
end

function solve_data()
   solve()
end

using Interpolations

#not working for now
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
    

    p = plot(Rvalues, Mvalues, legend = false, show = false)
    xlabel!(p, raw"Radius (km)")
    ylabel!(p, raw"Mass (M$_\odot$)")
    savefig("tov_plot.png")
end
