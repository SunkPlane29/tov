#NOTE: this is really a manual test, you should check the output files
#and see if they make sense and correspond to the description and comparisons
#with the source

using TOV
using Test
using Plots
using Printf

BASEPATH = ".."
outpath = joinpath(BASEPATH, "out")
eosfile = joinpath(BASEPATH, "test", "eos", "quarkmatter.csv")
eosheader = ["ρb", "p", "ϵ"]
eos = TOV.eos_from_file(eosfile, eosheader)

p₀ = 200

solution = TOV.solve_tov(p₀*TOV.MEVFM3_TO_PRESSURE_UNIT, eos, ϵsup=eos.eos_function(0.0))

plot(solution.r, solution.phi, xaxis=raw"$r$ (km)", yaxis=raw"$\phi$ (?)", label=false,
     title=@sprintf("p0 = %.2f MeV/fm^3", p₀))
savefig(joinpath(outpath, "deftest_rphi.png"))

plot(solution.r, solution.H, xaxis=raw"$r$ (km)", yaxis=raw"$H$ (?)", label=false,
     title=@sprintf("p0 = %.2f MeV/fm^3", p₀))
savefig(joinpath(outpath, "deftest_rH.png"))

plot(solution.r, solution.beta, xaxis=raw"$r$ (km)", yaxis=raw"$\beta$ (?)", label=false,
     title=@sprintf("p0 = %.2f MeV/fm^3", p₀))
savefig(joinpath(outpath, "deftest_rbeta.png"))

#TODO: figure out the units of these things I added
#TODO: also fix the problem with negative k2 and Lambda and also non-monotonic
#behavior in these two quantitities in function of mass

@printf("\nϵ_sup: %.4e MeV/fm^3\n", eos.eos_function(0.0)*TOV.PRESSURE_UNIT_TO_MEVFM3)
@printf("\nk2: %.4e \nLambda: %.4e\n", solution.k2, solution.Lambda)

p₀ = collect(range(1, 600, length=300)) .* TOV.MEVFM3_TO_PRESSURE_UNIT
sol = TOV.solve_sequence(p₀, eos, ϵsup=eos.eos_function(0.0), stepsize=1*TOV.SI_TO_LENGTH_UNIT)

plot(sol.M, sol.Lambda, xaxis=raw"$M$ (M$_{\odot}$)", yaxis=raw"$\Lambda$", label=false)
savefig(joinpath(outpath, "deftest_MLambda.png"))

plot(sol.M, sol.k2, xaxis=raw"$M$ (M$_{\odot}$)", yaxis=raw"$k_2$", label=false)
savefig(joinpath(outpath, "deftest_Mk2.png"))
