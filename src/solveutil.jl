#more generic function that calls solve_tov and catches errors that might appear, will also write the curve
#data to a .csv file (always named tov_data.csv NOTE: might change later) unless specified not to
#user must give an initial pressure and an equation of state (function of pressure) for the star
function solve(p₀::Real, eos::Function; write::Bool = true, stepsize::Real = 200*SI_TO_LENGTH_UNIT, n::Integer = 100000)::Curve
    curve = solve_tov(p₀, eos, stepsize, n)

    if write
        write_data(curve)
    end

    return curve
end

#FIXME: neutron_star giving totally innacurate results
#solve_plot calls solve and plots curves of pressure and mass in a single plot
function solve_plot(p₀::Real, eos::Function; stepsize::Real = 200*SI_TO_LENGTH_UNIT, n::Integer = 100000)
    curve = solve(p₀, eos, stepsize = stepsize, n = n)

    plot_curves(curve, raw"r(km)", raw"p(erg/cm³)", raw"M (M$_{\odot}$)", "single_star_plot.png")
end

#solve_data calls solve but does not plot, only writing the data to a file
function solve_data(p₀::Real, eos::Function; stepsize::Real = 200*SI_TO_LENGTH_UNIT, n::Integer = 100000)
    curve = solve(p₀, eos, stepsize = stepsize, n = n)

    return curve
end

#not working for now
function solve_star_curve(pa::Real, pb::Real, eos::Function; nstars::Integer = 1000, stepsize::Real = 200*SI_TO_LENGTH_UNIT, n::Integer = 100000)
    h = (pb - pa)/nstars

    Rvalues = []
    Mvalues = []

    for i = 1:nstars
        curve = solve(pa + (i-1)*h, eos, stepsize = stepsize, n = n, write = false)

        append!(Rvalues, last(curve.tvalues))
        append!(Mvalues, last(curve.yvalues))
    end

    p = plot(Rvalues, Mvalues, legend = false, show = false)
    xlabel!(p, raw"Radius (km)")
    ylabel!(p, raw"Mass (M$_\odot$)")
    savefig("tov_plot.png")
end
