using TOV, Test

@time begin
    @time @safetestset "Single Star Fermi Gas EoS" begin
        include("fermigas_single_test.jl")
    end
    @time @safetestset "Star Sequence Fermi Gas EoS" begin
        include("fermigas_sequence_test.jl")
    end
end
