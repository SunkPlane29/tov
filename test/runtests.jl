using Test
using DataFrames, CSV, Plots, Printf
using Revise
using TOV

#TODO: add Aqua

transform(col, val) = val
transform(col, val::Float64) = @sprintf("%.16e", val)

eos = EoS(joinpath(dirname(@__FILE__), "eos", "eos.csv"), ["ρb", "P", "ϵ"])
ϵ(P) = eos(P)

@testset "TOV" begin
    include("test_cublicspline.jl")
    include("test_diffsolve.jl")
    include("test_mrdiagram.jl")
    include("test_lovenumber.jl")
end