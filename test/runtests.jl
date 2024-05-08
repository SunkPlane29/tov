using Test
using CSV, DataFrames, DataInterpolations, Plots, Printf
using Revise
using TOV

#TODO: add Aqua

transform(col, val) = val
transform(col, val::Float64) = @sprintf("%.16e", val)
eosfile = joinpath(dirname(@__FILE__), "eos", "eos.csv")
df = CSV.File(eosfile, header=["ρb", "P", "ϵ"]) |> DataFrame
dfinterp = LinearInterpolation((df.ϵ).*MeVfm3, (df.P).*MeVfm3, extrapolate=true)
ϵ(P) = dfinterp(P)

@testset "TOV" begin
    include("test_diffsolve.jl")
    include("test_mrdiagram.jl")
    include("test_lovenumber.jl")
end