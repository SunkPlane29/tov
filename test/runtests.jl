using Test
using CSV, DataFrames
using Revise
using TOV

eosfile = joinpath(dirname(@__FILE__), "eos", "eos.csv")
df = CSV.File(eosfile, header=["ρb", "P", "ϵ"]) |> DataFrame

@testset "TOV" begin
    include("test_mrdiagram.jl")
    include("test_lovenumber.jl")
end