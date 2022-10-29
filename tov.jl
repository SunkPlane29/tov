include("diff.jl")
include("constants.jl")
include("eos/polytropic.jl")

using Printf

# pressure differential equation with some change for using
function pressure_equation(r::Real, p::Real, M::Real)::Real
    #TODO: maybe these checks should depend on debug mode or maybe they shouldn't even be here
    # if r == 0
    #     throw("r = 0 means division by 0 in pressure equation")
    # end

    # if M == 0
    #     throw("M = 0 means a 0 angular coefficient for the tangent, that can possibly cause the solving to stay at 0 throughout")
    # end

    @printf("radius:%.16e \npressure: %.16e \nmass: %.16e\n", r, p, M)

    #ATENTION: alpha should be defined in the code containing the eos
    return -α*ϵ_rel(p)*M/r^2
end

# is the actual mass coupled differential equation
function mass_equation(r::Real, p::Real, M::Real)::Real
    return β*(r^2)*ϵ_rel(p)
end

using QuadGK
# returns a mass value corresponding to the mass in the interior of the radius r
function mass_function(r::Real, p::Real, M::Real)::Real
    mm(r) = β*r^2*ϵ_rel(p)

    #TODO: not using err
    (int, err) = quadgk(mm, 0, r)
    return int
end

function solve_tov()::Tuple{Curve,Curve}
    p₀ = 1.0e-15
    r₀ = 1e-4
    m₀ = 1e-16
    stepsize = 100
    n = 5000

    return solve_system(pressure_equation, mass_equation, r₀, p₀, m₀, stepsize, n)
end
