using TOV, Test

@time begin
    @time @safetestset begin
        include("fermigas_single_test.jl")
    end
    @time @safetestset begin
        include("fermigas_sequence_test.jl")
    end
end
