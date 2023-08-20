#MANUAL TEST

include("../src/TOV.jl")
using .TOV
using Test

using Plots

# BASEPATH = dirname(dirname(pathof(TOV)))
BASEPATH = ".."
# eosfile = joinpath(BASEPATH, "test", "eos", "fermigas.csv")
# eosfile = joinpath(BASEPATH, "test", "eos", "fermigaslarge.csv")
eosfile = joinpath(BASEPATH, "test", "eos", "quarkmatter.csv")
# eosfile = joinpath(BASEPATH, "test", "eos", "cmf.csv")
# eosheader = ["p", "ϵ"]
eosheader = ["ρb", "p", "ϵ"]
# eosheader = ["T", "n_b", "Y_q", "p", "ϵ"]
eos = TOV.eos_from_file(eosfile, eosheader)

p₀ = eos.pressure
# p₀ = collect(range(eos.pressure[1], last(eos.pressure), length=300))

#1e-8 para fermigas
#1e-4 para quarkmatter or 5.82210288e-5 (maximum precision?)
#1e-7 para cmf
# sol = TOV.solve_sequence(p₀, eos, eps=1e-8)
sol = TOV.solve_sequence(p₀, eos, eps=1e-4)
# sol = TOV.solve_sequence(p₀, eos, eps=1e-7)

println(maximum(sol.M))

#NOTE: the error is most likelly in the solving of the TOV equations,
#specifically where I use a convergence cut-off
#NOTE: when settings eps for the cutoff, when eps ~ 1e-5 the mass starts
#to have artifacts just like the radius. Also, for more accurate maximum
#masses, we need more precise eps (of order ~ 1e-8)
#TODO: maybe I can get the order of magnitude of the pressure that
#has converged and use that same order of magnitude in all the other
#pressures (this would rely on channels I think)
#TODO: maybe I can let user choose between conversion in pressure or mass
#or even choose for a fixed stepsize method

plot(sol.p₀, sol.R, xlim=(-10, 1000), ylim=(0, 20), seriestype=:scatter,
     show=true, size=(1200, 900))
savefig("p0R.png")
plot(sol.p₀, sol.M, seriestype=:scatter, xlim=(-10, 1000), show=false,
     size=(1200, 900))
savefig("p0M.png")
plot(sol.R, sol.M, seriestype=:path, size=(1200, 900))
savefig("MR.png")
