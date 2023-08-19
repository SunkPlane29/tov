include("../src/TOV.jl")
using .TOV
using Test

using Plots

# BASEPATH = dirname(dirname(pathof(TOV)))
BASEPATH = ".."
eosfile = joinpath(BASEPATH, "test", "eos", "fermigas.csv")
eosheader = ["p", "ϵ"]
eos = TOV.eos_from_file(eosfile, eosheader)

p₀ = eos.pressure

sol = TOV.solve_sequence(p₀, eos, eps=1e-10)

M = []

for soli in sol
    println(last(soli.t) * TOV.LENGTH_UNIT_TO_SI * 1e-3)
    println(last(soli.u))
    append!(M, last(soli.u)[2])
end

plot(p₀ .* TOV.PRESSURE_UNIT_TO_MEVFM3, M)
