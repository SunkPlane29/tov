using TOV
using Test

#FIXME: NOT WORKING PROPERLY

BASEPATH = dirname(dirname(pathof(TOV)))
eosfile = joinpath(BASEPATH, "test", "eos", "fermigas.csv")
eosheader = ["p", "ϵ"]
eos = eos_from_file(eosfile, eosheader)

p₀ = eos.pressure

sol = solve_sequence(p₀, eos)
println(sol[100].t)
