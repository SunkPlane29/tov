# include("constants.jl")
# using .SIUnits

include("tov.jl")
include("eos/polytropic.jl")

include("util.jl")

#TODO: make solve function without saving data
function solve(p₀::Real, γ::Real, K::Real, write::Bool = true)::Curve
    # I think returning a one here would not change the equations (since its will not be added and only multiplied)
    eos(p) = p <= 0 ? 1 : polytrope(p, γ, K)

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

#suggestion use p₀ = 1.603e33 erg/cm³
function solve_plot(p₀::Real)
    #this make simpler to change from relativistic to non-relativistic later
    γ = γ_nonrel
    K = K_NONREL

    curve = solve(p₀, γ, K)

    plot_curves(curve, raw"r(km)", raw"p(dimensionless)", raw"M (M$_{\cdot}$)", "pressure_plot.png", "mass_plot.png")
end

function solve_data(p₀::Real)
    γ = γ_rel
    K = K_REL

    curve = solve(p₀, γ, K)

    return curve
end

#not working for now
function solve_star_curve(pa::Real, pb::Real)
    γ = γ_rel
    K = K_REL

    n = 10000
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
