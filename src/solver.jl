#I will define some units in which c = G = M⊙ = 1, i will also define convertion factors for
#the basics units in SI: m, kg, s
#
#1(unit of mass)    ≈ 1.989*10^30kg
#1(unit of length)  ≈ 1.477*10^3m
#1(unit of time)    ≈ 4.927*10^(-6)s

#TODO: refactor and add love number and potential
function solve_tov(p₀::Real, eos::Function, stepsize::Real ; n::Integer=100_000)::Curve
    r_init = 1e-8
    m_init = 1e-24
    p_init = p₀

    #equation of hydrostatic equilibrium
    pressure_eq(r, p, M) = begin
        newtonian = eos(p)*M/r^2
        special_rel_factor1 = 1 + p/eos(p)
        special_rel_factor2 = 1 + (4π*r^3*p)/M
        general_rel_factor = (1 - 2M/r)^(-1)
        slope = -newtonian*special_rel_factor1*special_rel_factor2*general_rel_factor
        r == 0 ? 0 : slope
    end
    #equation of mass continuity
    mass_eq(r, p, M) = begin
        slope = 4π*r^2*eos(p)
        r == 0 ? 0 : slope
    end
    #condition function required by the differential system solver
    condition_func(i, r, p, M) = p <= 0 || i > n ? false : true

    curve = Curve(Float64[], Float64[], Float64[])
    solve_system!(pressure_eq, mass_eq, r_init, p_init, m_init, curve, stepsize, condition_func)

    curve.tvalues = curve.tvalues*LENGTH_UNIT_TO_SI*1e-3
    curve.xvalues = curve.xvalues*PRESSURE_UNIT_TO_SI*JOULE_TO_MEV4*MEV4_TO_MEVFM3
    curve.yvalues = curve.yvalues
    return curve
end
