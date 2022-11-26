include("constants.jl")
using .SIUnits

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

#TODO: make solve function without saving data
function solve(p₀::Real, ϵ₀::Real, r₀::Real, γ::Real, K::Real, write::Bool = true)::Curve
    polytrp(p) = polytrope(p, γ, K)
    eos_const = ϵ₀_const(ϵ₀, γ)
    eos(p) = p <= 0 ? 0 : polytrp(p)*eos_const

    curve = try solve_tov(p₀, ϵ₀, r₀, eos)
        catch err
            println(err)
            return
        end

    if write
        write_data(curve)
    end

    return curve
end

#SUGESTION: use p₀ = 1e-16, ϵ₀ = 6e+38 and r₀ = 1e4 for relativistic limit (with p₀>1.54e-16)
#TODO: find parameters to fit non relativistic limit
#TODO: make some code later to pick these parameters to fit observations (? i don't know if it is allowed ?)
function solve_plot(p₀::Real, ϵ₀::Real, r₀::Real)
    #this make simpler to change from relativistic to non-relativistic later
    γ = γ_nonrel
    K = K_NONREL

    curve = solve(p₀, ϵ₀, r₀, γ, K)

    plot_functions(curve)
end

function solve_data(p₀::Real, ϵ₀::Real, r₀::Real)
    γ = γ_rel
    K = K_REL

    curve = solve(p₀, ϵ₀, r₀, γ, K)

    return curve
end

#not working for now
function solve_star_curve(pa::Real, pb::Real, ϵ₀::Real, r₀::Real)
    γ = γ_rel
    K = K_REL

    n = 10000
    h = (pb - pa)/n

    Rvalues = []
    Mvalues = []

    for i = 1:n
        curve = solve(pa + (i-1)*h, ϵ₀, r₀, γ, K)

        append!(Rvalues, last(curve.tvalues))
        append!(Mvalues, last(curve.yvalues))
    end
    

    p = plot(Rvalues, Mvalues, legend = false, show = false)
    xlabel!(p, raw"Radius (km)")
    ylabel!(p, raw"Mass (M$_\odot$)")
    savefig("tov_plot.png")
end
