#TODO: maybe put somewhere else
using Printf

#I will define some units in which c = G = M⊙ = 1, i will also define convertion factors for
#the basics units in SI: m, kg, s
#
#1(unit of mass)    = 1.989*10^30kg
#1(unit of length)  = 1.477*10^3m (also one half of the swartzchild radius of the sun)
#1(unit of time)    = 4.927*10^-6s

#solve_tov is a function that computes the solution of the Tollman-Oppenheimer-Volkoff (TOV) equations
#given a central pressure of the star and it's equation of state (EoS). stepsize and n (number of points)
#should also be passed as parameters since each star will have varying radiuses and masses
#the function returns a struct named curve that holds 3 vectors, one containing the array of radiuses, one
#containing the array of pressures and one containing the array of masses. The radius array should be in units
#of m or km (for now it's km), the pressure array can be in any units (but one could use units that help
#with comparison with other sources) and the mass array should always be in units of solar mass
function solve_tov(p₀::Real, eos::Function, stepsize::Real, n::Integer)::Curve
    #TODO: i had to choose these initial values very carefully, also, they might affect the solution a bit
    #too much
    #NOTE: from what I quickly analised it does change the solution, but the change is not that significant,
    #the change will be like an error
    r_init = 1e-8
    m_init = 1e-24
    p_init = p₀

    #TODO: print statements should not be on modules and non main files
    @printf("Solving TOV with p₀ = %.8e\n", p₀)

    #equation of hydrostatic equilibrium with exception thrown when the pressure is negative
    pressure_eq(r, p, M) = begin
        #TODO: think of better error management
        if p <= 0
            throw(ErrorException("error: pressure is lesser or equal to zero"))
        end

        newtonian = eos(p)*M/r^2
        special_rel_factor1 = 1 + p/eos(p)
        special_rel_factor2 = 1 + (4π*r^3*p)/M
        general_rel_factor = (1 - 2M/r)^(-1)
        slope = -newtonian#*special_rel_factor1*special_rel_factor2*general_rel_factor
        r == 0 ? 0 : slope
    end
    #equation of mass continuity
    mass_eq(r, p, M) = begin
        if p <= 0
            throw(ErrorException("error: pressure is lesser or equal to zero"))
        end

        slope = 4π*r^2*eos(p)
        r == 0 ? 0 : slope
    end
    #condition function required by the differential system solver
    condition_func(i, r, p, M) = p <= 0 || i > 100000 ? false : true

    curve = Curve(Float64[], Float64[], Float64[])
    try
        solve_system!(pressure_eq, mass_eq, r_init, p_init, m_init, curve, stepsize, condition_func)
    catch e
        @printf("error while solving diff equations: %s\n", e)
        curve.tvalues = curve.tvalues*LENGTH_UNIT_TO_SI*1e-3
        curve.xvalues = curve.xvalues*PRESSURE_UNIT_TO_SI*1e1#*SI_TO_GEV_FM3
        curve.yvalues = curve.yvalues
        return curve
    end

    @printf("expected error ;)\n")
    return curve
end
