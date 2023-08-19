include("../src/TOV.jl")
using .TOV
# using TOV
using Test

struct ExpectedSolution
    p₀::Real
    p::Real
    R::Real
    M::Real
end

# For p₀=1.0MeV/fm³ we have R=19.884000014770535km and M=0.3295759365485324MSOLAR
# For p₀=10.0MeV/fm³ we have R=14.869000014766742km and M=0.5344830305160919MSOLAR
# For p₀=25.0MeV/fm³ we have R=13.054000014765366km and M=0.6133481871959284MSOLAR

# BASEPATH = dirname(dirname(pathof(TOV)))
BASEPATH = ".."
eosfile = joinpath(BASEPATH, "test", "eos", "fermigas.csv")
eosheader = ["p", "ϵ"]
eos = TOV.eos_from_file(eosfile, eosheader)

want = ExpectedSolution(1.0, 0.0, 19.884000014770535, 0.3295759365485324)
solution = TOV.solve_tov(want.p₀*TOV.MEVFM3_TO_PRESSURE_UNIT, eos, eps=1e-10)

p = last(solution.p)
R = last(solution.r)
M = last(solution.M)

@test want.p ≈ p atol=1e-2
#NOTE: this tolerance may be due to the previous method of solving the TOV equation
#being already innacurate, so this method might actually be better
@test want.R ≈ R atol=4e-1
@test want.M ≈ M atol=1e-2

want = ExpectedSolution(10.0, 0.0, 14.869000014766742, 0.5344830305160919)
solution = TOV.solve_tov(want.p₀*TOV.MEVFM3_TO_PRESSURE_UNIT, eos, eps=1e-10)

p = last(solution.p)
R = last(solution.r)
M = last(solution.M)

@test want.p ≈ p atol=1e-2
@test want.R ≈ R atol=4e-1
@test want.M ≈ M atol=1e-2

want = ExpectedSolution(25.0, 0.0, 13.054000014765366, 0.6133481871959284)
solution = TOV.solve_tov(want.p₀*TOV.MEVFM3_TO_PRESSURE_UNIT, eos, eps=1e-9)

p = last(solution.p)
R = last(solution.r)
M = last(solution.M)

@test want.p ≈ p atol=1e-2
    @test want.R ≈ R atol=4e-1
@test want.M ≈ M atol=1e-2
