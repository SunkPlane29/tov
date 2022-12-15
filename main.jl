include("constants.jl")

include("tov.jl")
include("eos/polytropic_whitedwarf.jl")

include("util.jl")

function solve(p₀::Real, γ::Real, K::Real, write::Bool = true)::Curve
    eos(p) = polytrope(p, γ, K)

    curve = try solve_tov(p₀, eos)
        catch err
            println(err)
            return
        end

    if write
        write_data(curve)
    end

    return curve
end

#TODO: i don't know if there is something wrong with these units or the equation of state, but i get different
#results from the paper (apparently only on the relativistic neutron star verion, might be a problem with the eos)
#suggestion use p₀ = 7.463e22 Pa J/m³ for relativistic version white dwarf (and higher power)
#suggestion use p₀ = 2.488e23 J/m³ for non-relativistic version white dwarf (and higher power)
#FIXME: non-relativistic white dwarf giving totally different results for the pressure given
function solve_plot(p₀::Real, relativistic::Bool = true)
    γ = relativistic ? γ_rel : γ_nonrel
    K = relativistic ? K_REL : K_NONREL

    curve = solve(p₀, γ, K)

    plot_curves(curve, raw"r(km)", raw"p(MeV/fm³)", raw"M (M$_{\cdot}$)", "pressure_plot.png", "mass_plot.png")
end

function solve_data(p₀::Real)
    γ = relativistic ? γ_rel : γ_nonrel
    K = relativistic ? K_REL : K_NONREL

    curve = solve(p₀, γ, K)

    return curve
end

#not working for now
function solve_star_curve(pa::Real, pb::Real)
    γ = γ_nonrel
    K = K_NONREL

    n = 1000
    h = (pb - pa)/n

    Rvalues = []
    Mvalues = []

    for i = 1:n
        curve = solve(pa + (i-1)*h, γ, K)

        append!(Rvalues, last(curve.tvalues))
        append!(Mvalues, last(curve.yvalues))
    end
    

    p = plot(Rvalues, Mvalues, legend = false, show = false)
    xlabel!(p, raw"Radius (km)")
    ylabel!(p, raw"Mass (M$_\odot$)")
    savefig("tov_plot.png")
end
