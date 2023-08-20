#MANUAL TEST

include("../src/TOV.jl")
using .TOV
using Test

using Plots

# BASEPATH = dirname(dirname(pathof(TOV)))
BASEPATH = ".."
eosfile = joinpath(BASEPATH, "test", "eos", "fermigas.csv")
# eosfile = joinpath(BASEPATH, "test", "eos", "fermigaslarge.csv")
# eosfile = joinpath(BASEPATH, "test", "eos", "quarkmatter.csv")
# eosfile = joinpath(BASEPATH, "test", "eos", "cmf.csv")
eosheader = ["p", "ϵ"]
# eosheader = ["ρb", "p", "ϵ"]
# eosheader = ["T", "n_b", "Y_q", "p", "ϵ"]
eos = TOV.eos_from_file(eosfile, eosheader)

p₀ = eos.pressure

#1e-7 para fermigas
#6e-5 para quarkmatter
#1e-7 para cmf
# sol = TOV.solve_sequence(p₀, eos, eps=6e-5)
sol = TOV.solve_sequence(p₀, eos, eps=1e-7)

#TODO: maybe I can get the order of magnitude of the pressure that
#has converged and use that same order of magnitude in all the other
#pressures (this would rely on channels I think)

plot(sol.p₀, sol.R, xlim=(-10, 250), ylim=(5, 20), seriestype=:scatter,
     show=true, size=(1200, 900))
savefig("p0R.png")
plot(sol.p₀, sol.M, seriestype=:scatter, xlim=(-10, 1000), show=false,
     size=(1200, 900))
savefig("p0M.png")
plot(sol.R, sol.M, seriestype=:path, size=(1200, 900))
savefig("MR.png")
