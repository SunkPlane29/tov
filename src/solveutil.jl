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

using CSV
using DataFrames

#not working for now
#NOTE: this function uses all threads avaliable in julia program (a library shouldn't do multithreaded code, but
#I think this solution is suitable here, since it reduces in 50%+ the execution time)
#TODO: there might be a solution using GPU, but this solution can be in another function
#solve_mrdiagram solves TOV equations for a range of initial pressures. the resulting plot is a Mass Radius diagram for
#the chosen star type. Note that this function uses multithreaded code when available, i.e., when the julia repl or inter
#preter is called with --threads (nthreads) argument
function solve_mrdiagram(pa::Real, pb::Real, eos::Function; nstars::Integer = 1000, stepsize::Real = 200*SI_TO_LENGTH_UNIT, n::Integer = 100000)::Curve
    h = (pb - pa)/nstars

    pvalues = []
    Rvalues = []
    Mvalues = []

    #multithreaded for loop
    l = ReentrantLock()
    Threads.@threads :dynamic for i = 1:nstars
        p₀ = pa + (i-1)*h
        curve = solve(p₀, eos, stepsize = stepsize, n = n, write = false)

        lock(l)
        try
            append!(pvalues, p₀)
            append!(Rvalues, last(curve.tvalues))
            append!(Mvalues, last(curve.yvalues))
        finally
            unlock(l)
        end
    end

    #since the tov was solved in a multithreaded way, the resulting M and R values are unordered, to solve this
    #one can simply apply a sort in p₀ values and then apply the same sorting on M and R values
    perm = sortperm(pvalues, alg = QuickSort)
    pvalues = pvalues[perm]
    Rvalues = Rvalues[perm]
    Mvalues = Mvalues[perm]

    df = DataFrame()
    df.radius = Rvalues
    df.mass = Mvalues
    CSV.write("mrdiagram.csv", df)

    p = plot(Rvalues, Mvalues, legend = false, show = false)
    xlabel!(p, raw"Radius (km)")
    ylabel!(p, raw"Mass (M$_\odot$)")
    savefig("mrdiagram.png")

    return Curve(pvalues, Rvalues, Mvalues)
end
