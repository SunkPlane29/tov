#NOTE: this is really a manual test, you should check the output files
#and see if they make sense and correspond to the description and comparisons
#with the source

using TOV
using Test
using Plots
using Printf
using DataFrames
using CSV

BASEPATH = ".."
outpath = joinpath(BASEPATH, "out")
eosfile1 = joinpath(BASEPATH, "test", "eos", "quarkmatter.csv")
eosfile2 = joinpath(BASEPATH, "test", "eos", "cmf.csv")
eosfile3 = joinpath(BASEPATH, "test", "eos", "su2njl_eos.csv")
eosheader1 = ["ρb", "p", "ϵ"]
eosheader2 = ["T", "n_b", "Y_q", "p", "ϵ"]
eosheader3 = ["ϵ", "p", "n_b", "μ_b"]
eos1 = TOV.eos_from_file(eosfile1, eosheader1)
eos2 = TOV.eos_from_file(eosfile2, eosheader2)
eos3 = TOV.eos_from_file(eosfile3, eosheader3)

p₀ = 200

eos = eos1

ϵsup = eos.eos_function(0.0)

sol = TOV.solve_tov(p₀*TOV.MEVFM3_TO_PRESSURE_UNIT, eos, ϵsup=ϵsup)

plot(sol.r, sol.phi, xaxis=raw"$r$ (km)", yaxis=raw"$\phi$ (?)", label=false,
     title=@sprintf("p0 = %.2f MeV/fm^3", p₀))
savefig(joinpath(outpath, "deftest_rphi.png"))

plot(sol.r, sol.H, xaxis=raw"$r$ (km)", yaxis=raw"$H$ (?)", label=false,
     title=@sprintf("p0 = %.2f MeV/fm^3", p₀))
plot!(sol.r, sol.r .^ 2, label="Initial condition", linestyle=:dash)
savefig(joinpath(outpath, "deftest_rH.png"))

plot(sol.r, sol.beta, xaxis=raw"$r$ (km)", yaxis=raw"$\beta$ (?)", label=false,
     title=@sprintf("p0 = %.2f MeV/fm^3", p₀))
plot!(sol.r, 2 .* sol.r, label="Initial condition", linestyle=:dash)
savefig(joinpath(outpath, "deftest_rbeta.png"))

#TODO: also fix the problem with negative k2 and Lambda and also non-monotonic
#behavior in these two quantitities in function of mass

@printf("\nϵ_sup: %.4e MeV/fm^3\n", eos.eos_function(0.0)*TOV.PRESSURE_UNIT_TO_MEVFM3)
@printf("\nk2: %.4e \nLambda: %.4e\n", sol.k2, sol.Lambda)

# p₀ = collect(range(1, 600, length=300)) .* TOV.MEVFM3_TO_PRESSURE_UNIT
p₀ = [collect(range(1, 6, length=50)); collect(range(6, 2000, length=300))] .* TOV.MEVFM3_TO_PRESSURE_UNIT
sol = TOV.solve_sequence(p₀, eos, ϵsup=ϵsup, stepsize=1*TOV.SI_TO_LENGTH_UNIT)

plot(sol.R, sol.M, xaxis=raw"$R$ (km)", yaxis=raw"$M$ (M$_{\odot}$)", label=false)
savefig(joinpath(outpath, "deftest_MR.png"))

plot(sol.M, sol.Lambda, xaxis=raw"$M$ (M$_{\odot}$)", yaxis=raw"$\log \Lambda$", label=false, 
yscale=:log10, yticks = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10])
savefig(joinpath(outpath, "deftest_MLambda.png"))

plot(sol.M, sol.k2, xaxis=raw"$M$ (M$_{\odot}$)", yaxis=raw"$k_2$", label=false)
savefig(joinpath(outpath, "deftest_Mk2.png"))

plot(sol.R, sol.H, xaxis=raw"$R$ (km)", yaxis=raw"$H$ (?)", label=false)
savefig(joinpath(outpath, "deftest_RH.png"))

plot(sol.R, sol.Beta, xaxis=raw"$R$ (km)", yaxis=raw"$\Beta$ (?)", label=false)
savefig(joinpath(outpath, "deftest_RBeta.png"))

plot(sol.M, sol.H, xaxis=raw"$M$ (M$_{\odot}$)", yaxis=raw"$H$ (?)", label=false)
savefig(joinpath(outpath, "deftest_MH.png"))

plot(sol.M, sol.Beta, xaxis=raw"$M$ (M$_{\odot}$)", yaxis=raw"$\Beta$ (?)", label=false)
savefig(joinpath(outpath, "deftest_MBeta.png"))

mr = DataFrame(p0=sol.p₀, R=sol.R, M=sol.M, k2=sol.k2, Lambda=sol.Lambda, Phi=sol.Phi,
               H=sol.H, Beta=sol.Beta)
CSV.write(joinpath(outpath, "deftest_MR.csv"), mr)
