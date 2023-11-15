#MANUAL TEST, JUST CHECK THE RESULTS
#
#MAXIMUM MASS OF THE FIRST EOS SHOULD BE AROUND 1.97 SOLARMASS
#MAXIMUM MASS OF THE SECOND EOS SHOULD BE AROUND 2.02 SOLARMASS

using TOV
using Test

using Plots
using Printf

using CSV
using DataFrames

BASEPATH = dirname(dirname(pathof(TOV)))
run(`mkdir out -p`)
out1 = joinpath(BASEPATH, "out", "p0MR1.csv")
out2 = joinpath(BASEPATH, "out", "p0MR2.csv")
out3 = joinpath(BASEPATH, "out", "MR1.png")
out4 = joinpath(BASEPATH, "out", "MR2.png")
run(`touch $out1 $out2 $out3 $out4`)

eosfile1 = joinpath(BASEPATH, "test", "eos", "quarkmatter.csv")
eosfile2 = joinpath(BASEPATH, "test", "eos", "cmf.csv")
eosheader1 = ["ρb", "p", "ϵ"]
eosheader2 = ["T", "n_b", "Y_q", "p", "ϵ"]

@time begin
    eos1 = TOV.eos_from_file(eosfile1, eosheader1)
    p01 = collect(range(1, 600, length=300)) .* TOV.MEVFM3_TO_PRESSURE_UNIT

    sol = TOV.solve_sequence(p01, eos1, stepsize=1*TOV.SI_TO_LENGTH_UNIT)

    CSV.write(joinpath(BASEPATH, "out", "seqtest_p0MR1.csv"), DataFrame(p0=sol.p₀, M=sol.M, R=sol.R))

    plot(sol.R, sol.M, seriestype=:path,
         title=@sprintf("Maximum mass: %.4f", maximum(sol.M)), label=false,
         yaxis=raw"$M$ (M$_{\odot}$)", xaxis=raw"$R$ (km)")

    savefig(joinpath(BASEPATH, "out", "seqtest_MR1.png"))

    @test 1.9722 ≈ maximum(sol.M) atol=1e-2
end

@time begin
    eos2 = TOV.eos_from_file(eosfile2, eosheader2)
    p02 = [collect(range(1, 6, length=50)); collect(range(6, 600, length=250))] .* TOV.MEVFM3_TO_PRESSURE_UNIT

    sol = TOV.solve_sequence(p02, eos2, stepsize=1*TOV.SI_TO_LENGTH_UNIT)

    CSV.write(joinpath(BASEPATH, "out", "seqtest_p0MR2.csv"), DataFrame(p0=sol.p₀, M=sol.M, R=sol.R))

    plot(sol.R, sol.M, seriestype=:path,
         title=@sprintf("Maximum mass: %.4f", maximum(sol.M)), label=false,
         yaxis=raw"$M$ (M$_{\odot}$)", xaxis=raw"$R$ (km)")

    savefig(joinpath(BASEPATH, "out", "seqtest_MR2.png"))

    @test 2.0203 ≈ maximum(sol.M) atol=1e-2
end
