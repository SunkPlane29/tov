#MANUAL TEST

include("../src/TOV.jl")
using .TOV
using Test

using Plots

# BASEPATH = dirname(dirname(pathof(TOV)))
BASEPATH = ".."
# eosfile = joinpath(BASEPATH, "test", "eos", "fermigas.csv")
# eosfile = joinpath(BASEPATH, "test", "eos", "fermigaslarge.csv")
# eosfile = joinpath(BASEPATH, "test", "eos", "quarkmatter.csv")
eosfile = joinpath(BASEPATH, "test", "eos", "cmf.csv")
# eosheader = ["p", "ϵ"]
# eosheader = ["ρb", "p", "ϵ"]
# eosheader = ["T", "n_b", "Y_q", "p", "ϵ"]
eos = TOV.eos_from_file(eosfile, eosheader)

p₀ = eos.pressure

#1e-7 para fermigas
#6e-5 para quarkmatter
sol = TOV.solve_sequence(p₀, eos, eps=1e-7)

#TODO: as is it solves MR diagrams pretty well, besides the problem
#of having the plot filled with artifacts

println(sol.R)
println(sol.M)

# plot(sol.p₀, sol.R, xlim=(-10, 500), show=true)
# plot(sol.p₀, sol.M, xlim=(-10, 1000), show=false)
plot(sol.R, sol.M)
