using Test
using CSV, DataFrames
using Revise
using TOV

@testset "TOV" begin
    include("test_mrdiagram.jl")
    include("test_lovenumber.jl")
end