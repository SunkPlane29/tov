@testset "MR Diagram" begin

P0 = range(1MeVfm3, 600MeVfm3, length=200)

mrdiagram = solvemrdiagram(P0, ϵ)
plot(mrdiagram[:, 3], mrdiagram[:, 2], label="TOV", ylabel="M [M⊙]", xlabel="R [km]", legend=:topleft)
savefig(joinpath(dirname(@__FILE__), "out", "test_mrdiagram.png"))
outdf = DataFrame(P0=P0, M=mrdiagram[:, 2], R=mrdiagram[:, 3])
CSV.write(joinpath(dirname(@__FILE__), "out", "test_mrdiagram.csv"), outdf)

end