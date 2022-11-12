include("diff.jl")
include("constants.jl")
include("eos/polytropic.jl")

using Printf

function solve_tov(p_init::Real)::Curve
    r_init = 0
    m_init = 0
    #TODO: I think that if i ajust very VERY well this one I may not have to change it never again
    r₀ = 0.9319e4 #with p_init = 1.54e-16 and r₀ = 0.319e4 we get very near the results in the paper (varying the exponent in p_init also get near the other results in the paper)
    n = 10000
    stepsize = (r₀ - r_init)/n

    @printf("Solving TOV with p_init = %.8e\n", p_init)

    R₀ = G*MSOLAR/c^2
    α = R₀/r₀
    ϵ₀ = ((1/K_REL)*(R₀/(α*r₀))^γ_rel)^(1/(γ_rel-1))
    β = (4π*r₀^3*ϵ₀)/(MSOLAR*c^2)

    @printf("Constants α = %.8e, β = %.8e, ϵ₀ = %.8e\n", α, β, ϵ₀)

    eos(p) = p < 0 ? 0 : p^(1/γ_rel)

    #TODO: think later about refactoring this
    #TODO: also, too much clustered code, but this equations need to get some constants from environtment
    pressure_eq(r, p, M) = r == 0 ? 0 : -α*((M*eos(p))/r^2)
    mass_eq(r, p, M) = r == 0 ? 0 : β*r^2*eos(p)
    condition_func(i, r, p, M) = p <= 0 || i > 100000 ? false : true

    curve = solve_system(pressure_eq, mass_eq, r_init, p_init, m_init, stepsize, condition_func)
    #TODO: check these convertions later
    curve.tvalues = curve.tvalues*r₀*1e-3
    curve.xvalues = curve.xvalues*ϵ₀
    return curve
end
