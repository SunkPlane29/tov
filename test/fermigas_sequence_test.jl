#MANUAL TEST

include("../src/TOV.jl")
using .TOV
using Test

using Plots
using Printf

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

#NOTE: using the same p₀ (without interpolation) leads probably to
#lesser errors
p₀ = eos.pressure
# p₀ = collect(range(eos.pressure[1], last(eos.pressure), length=300))

singlesol = TOV.solve_tov(200*TOV.MEVFM3_TO_PRESSURE_UNIT, eos, eps=1e-8)
plot(singlesol.r, singlesol.p,
     title=@sprintf("R: %.4f, M: %.4f, p₀: %.4e", last(singlesol.r),
                    last(singlesol.M), last(singlesol.p₀)))
savefig("rp.png")
plot(singlesol.r, singlesol.M)
savefig("rM.png")

#1e-8 para fermigas (1e-7 for large)
#1e-4 para quarkmatter or 5.82210288e-5 (maximum precision?)
#1e-7 para cmf
# sol = TOV.solve_sequence(p₀, eos, eps=1e-8)
# sol = TOV.solve_sequence(p₀, eos, eps=1e-7)
sol = TOV.solve_sequence(p₀, eos, eps=1e-4)
# sol = TOV.solve_sequence(p₀, eos, eps=1e-7)

println(maximum(sol.R))
println(maximum(sol.M))

#NOTE: the error is most likelly in the solving of the TOV equations,
#specifically where I use a convergence cut-off
#NOTE: when settings eps for the cutoff, when eps ~ 1e-5 the mass starts
#to have artifacts just like the radius. Also, for more accurate maximum
#NOTE: I was able to set a minumum to the adaptive stepsize of the RK
#algorithm but the numeric errors are still there, why?
#masses, we need more precise eps (of order ~ 1e-8)
#NOTE: Minimum stepsize helped in the case of quarkmatter eos but made it
#worse in the case of cmf eos and also in the case of fermigas eos
#NOTE: stepsize in quarkmatter eos might just be off, because the equation
#is so stiff, it might be overdoing the step
#TODO: the error might be in interpolation?

plot(sol.p₀, sol.R, xlim=(-10, 1000), ylim=(0, 20), seriestype=:scatter,
     show=true, size=(1200, 900))
savefig("p0R.png")
plot(sol.p₀, sol.M, seriestype=:scatter, xlim=(-10, 1000), show=false,
     size=(1200, 900))
savefig("p0M.png")
plot(sol.R, sol.M, seriestype=:path, size=(1200, 900))
savefig("MR.png")
