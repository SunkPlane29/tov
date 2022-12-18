function solve(p₀::Real, eos::Function; write::Bool = true, stepsize::Real = 200*SI_TO_LENGTH_UNIT, n::Integer = 100000)::Curve
    curve = try solve_tov(p₀, eos, stepsize, n)
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
#suggestion use p₀ = 7.463e22 J/m³ for relativistic version white dwarf (and higher power)
#suggestion use p₀ = 2.488e23 J/m³ for non-relativistic version white dwarf (and higher power)
#FIXME: non-relativistic white dwarf giving totally different results for the pressure given
function solve_plot(p₀::Real, eos::Function)
    curve = solve(p₀, eos)

    plot_curves(curve, raw"r(km)", raw"p(MeV/fm³)", raw"M (M$_{\cdot}$)", "pressure_plot.png", "mass_plot.png")
end

function solve_data(p₀::Real, eos::Function)
    curve = solve(p₀, eos)

    return curve
end

#not working for now
function solve_star_curve(pa::Real, pb::Real, eos::Function)
    n = 1000
    h = (pb - pa)/n

    Rvalues = []
    Mvalues = []

    for i = 1:n
        curve = solve(pa + (i-1)*h, eos, write = false)

        append!(Rvalues, last(curve.tvalues))
        append!(Mvalues, last(curve.yvalues))
    end
    

    p = plot(Rvalues, Mvalues, legend = false, show = false)
    xlabel!(p, raw"Radius (km)")
    ylabel!(p, raw"Mass (M$_\odot$)")
    savefig("tov_plot.png")
end
