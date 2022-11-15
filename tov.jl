include("diff.jl")

using Printf

#FIXME: or maybe not. But in julia there is a problem of using modules and defining them
#so for this function to work the constants need to be included from constants.jl and also
#the system of units need to be declared with "using .$System$Units"
function solve_tov(p₀::Real, ϵ₀::Real, r₀::Real, eos::Function)::Curve
    r_init = 0
    m_init = 0
    p_init = p₀
    n = 10000
    stepsize = (r₀ - r_init)/n

    @printf("Solving TOV with p_init = %.8e\n", p₀)

    R₀ = G*MSOLAR/c^2
    α = R₀/r₀
    β = (4π*r₀^3*ϵ₀)/(MSOLAR*c^2)

    @printf("Constants α = %.8e, β = %.8e", α, β)

    pressure_eq(r, p, M) = r == 0 ? 0 : -α*((M*eos(p))/r^2)
    mass_eq(r, p, M) = r == 0 ? 0 : β*r^2*eos(p)
    condition_func(i, r, p, M) = p <= 0 || i > 100000 ? false : true

    curve = solve_system(pressure_eq, mass_eq, r_init, p_init, m_init, stepsize, condition_func)

    #convertion from dimensionless ̄r and ̄p and also expressing radius in km
    #TODO: the convertion to km is not flexible in terms of system of units
    curve.tvalues = curve.tvalues*r₀*1e-3
    curve.xvalues = curve.xvalues*ϵ₀
    return curve
end
