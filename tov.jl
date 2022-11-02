include("diff.jl")
include("constants.jl")
include("eos/polytropic.jl")

using Printf

α = G*MSOLAR/c^2
β = 4π/c^2
#β = 4π/MSOLAR*c^2

# pressure differential equation with some change for using
function pressure_equation(r::Real, p::Real, M::Real)::Real
    #TODO: maybe these checks should depend on debug mode or maybe they shouldn't even be here
    # if r == 0
    #     throw("r = 0 means division by 0 in pressure equation")
    # end

    # if M == 0
    #     throw("M = 0 means a 0 angular coefficient for the tangent, that can possibly cause the solving to stay at 0 throughout")
    # end

    @printf("r: %.8e, p: %.8e, ϵ_rel(): %.8e, M: %.8e\n", r, p, ϵ_rel(p), M)
    println(-α*ϵ_rel(p)*mass_function(r, p, M)/r^2)
    #TODO: this return value is possibly a very small number and that causes p to not shrink enougth)
    #FIXME: ϵ_rel < p. WHYYYY???
    #return (-α*ϵ_rel(p)*M/r^2)#*(1 + p/ϵ_rel(p))*(1 + (4π*r^3*p/mass_function(r, p, M)*c^2))*((1 - (2*G*mass_function(r, p, M)/c^2*r))^(-1))
    return (-α*ϵ_rel(p)*mass_function(r, p, M)/r^2)#*(1 + p/ϵ_rel(p))*(1 + (4π*r^3*p/mass_function(r, p, M)*c^2))*((1 - (2*G*mass_function(r, p, M)/c^2*r))^(-1))
end

# is the actual mass coupled differential equation
function mass_equation(r::Real, p::Real, M::Real)::Real
    return β*(r^2)*ϵ_rel(p)
end

using QuadGK
# returns a mass value corresponding to the mass in the interior of the radius r
function mass_function(r::Real, p::Real, M::Real)::Real
    mm(r) = β*r^2*ϵ_rel(p) ;

    #TODO: not using err
    (int, err) = quadgk(mm, 0, r)
    return int
end

function solve_tov(p₀::Real)::Curve
    r₀ = 1.0e-16
    m₀ = 1.0e-16
    rmax = 10000
    n = 100
    stepsize = (rmax - r₀)/n

    @printf("Solving TOV with p₀ = %.8e\n", p₀)

    curve = solve_system(pressure_equation, mass_equation, r₀, p₀, m₀, stepsize, n)
    return curve
end
