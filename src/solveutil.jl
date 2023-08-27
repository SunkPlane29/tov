#NOTE: these utility functions write the results into datafiles and plot images, a library code shouldn't
#do these sort of things (also use multithreaded), but..., I think this solution is suitable here since
#this wasn't made to be run on servers or in critical applications
#
#more generic function that calls solve_tov and catches errors that might appear, will also write the curve
#data to a .csv file unless specified not to
#user must give an initial pressure and an equation of state (function of pressure) for the star
function solve(p₀::Real, eos::Function;
               write::Bool=true, stepsize::Real=1*SI_TO_LENGTH_UNIT, n::Integer=100_000)::Curve
    curve = solve_tov(p₀, eos, stepsize, n=n)

    if write
        write_data(curve)
    end

    return curve
end

#solve_plot calls solve and plots curves of pressure and mass in a single plot
function solve_plot(p₀::Real, eos::Function;
                    stepsize::Real=1*SI_TO_LENGTH_UNIT, n::Integer=100_000, name::String="single_star_plot.png")
    curve = solve(p₀, eos, stepsize=stepsize, n=n)::Curve

    plot_curves(curve, raw"r (km)", raw"p (MeV/fm³)", raw"M (M$_{\odot}$)", name)
end

#solve_data calls solve but does not plot, only writing the data to a file
function solve_data(p₀::Real, eos::Function;
                    stepsize::Real=1*SI_TO_LENGTH_UNIT, n::Integer=100_000)::Curve
    curve = solve(p₀, eos, stepsize=stepsize, n=n)

    return curve
end

using CSV
using DataFrames

#solve_mrdiagram solves TOV equations for a range of initial pressures. the resulting plot is a Mass Radius diagram for
#the chosen star type. Note that this function uses multithreaded code when available, i.e., when the julia repl or inter
#preter is called with --threads (nthreads) argument
#NOTE: this function uses all threads avaliable in julia program (a library shouldn't do multithreaded code, but
#I think this solution is suitable here, since it reduces in 50%+ the execution time)
function solve_mrdiagram(pa::Real, pb::Real, eos::Function;
                         write_csv::Bool=true, plot::Bool=true, csvname::String="mrdiagram.csv",
                         plotname::String="mrdiagram.csv", nstars::Integer=500, stepsize::Real=1*SI_TO_LENGTH_UNIT,
                         n::Integer=100_000)::Curve

    h = (pb - pa)/nstars

    pvalues = []
    Rvalues = []
    Mvalues = []

    #multithreaded for loop
    l = ReentrantLock()
    Threads.@threads :dynamic for i = 1:nstars
        p₀ = pa + (i-1)*h
        curve = solve(p₀, eos, stepsize=stepsize, write=false)

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
    perm = sortperm(pvalues, alg=QuickSort)
    pvalues = pvalues[perm]
    Rvalues = Rvalues[perm]
    Mvalues = Mvalues[perm]

    if write_csv
        df = DataFrame()
        df.p0 = pvalues
        df.radius = Rvalues
        df.mass = Mvalues
        CSV.write(csvname, df)
    end

    if plot
        p = plot(Rvalues, Mvalues, legend=false, show=false, xaxis=raw"R (km)", yaxis=raw"M (M$_\odot$)")
            savefig(p, plotname)
    end

    return Curve(pvalues, Rvalues, Mvalues)
end
