using Test
using CSV, DataFrames, DataInterpolations, Plots
using Revise
using TOV

eosfile = joinpath(dirname(@__FILE__), "eos", "eos.csv")
df = CSV.File(eosfile, header=["ρb", "P", "ϵ"]) |> DataFrame
dfinterp = LinearInterpolation((df.ϵ).*MeVfm3, (df.P).*MeVfm3, extrapolate=true)
ϵ(P) = dfinterp(P)

@testset "TOV" begin
    include("test_diffsolve.jl")
    include("test_mrdiagram.jl")
    include("test_lovenumber.jl")
end