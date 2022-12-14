include("diff.jl")
include("constants.jl")

using Printf

#I will define some units in which c = G = M⊙ = 1, i will also define convertion factors for
#the basics units in SI: m, kg, s
#
#1(unit of mass)    = 1.989*10^30kg
#1(unit of length)  = 1.477*10^3m (also one half of the swartzchild radius of the sun)
#1(unit of time)    = 4.927*10^-6s

using Debugger

function solve_tov(p₀::Real, eos::Function)::Curve
    #TODO: i had to choose these initial values very carefully, also, they might affect the solution a bit
    #too much
    #NOTE: from what I quickly analised it does change the solution, but the change is not that significant,
    #the change will be like an error
    r_init = 1e-8
    m_init = 1e-24
    p_init = p₀
    n = 10000
    stepsize = 150*SI_TO_LENGTH_UNIT#TODO: change later to argument

    @printf("Solving TOV with p₀ = %.8e\n", p₀)

    pressure_eq(r, p, M) = begin
        if p <= 0
            throw(ErrorException("error: pressure is lesser or equal to zero"))
        end

        newtonian = eos(p)*M/r^2
        special_rel_factor1 = 1 + p/eos(p)
        special_rel_factor2 = 1 + (4π*r^3*p)/M
        general_rel_factor = (1 - 2M/r)^(-1)
        slope = -newtonian#*special_rel_factor1*special_rel_factor2*general_rel_factor
        @bp
        r == 0 ? 0 : slope
    end
    mass_eq(r, p, M) = begin
        if p <= 0
            throw(ErrorException("error: pressure is lesser or equal to zero"))
        end

        slope = 4π*r^2*eos(p)
        r == 0 ? 0 : slope
    end
    condition_func(i, r, p, M) = p <= 0 || i > 100000 ? false : true

    curve = Curve(Real[], Real[], Real[])

    try
        solve_system!(pressure_eq, mass_eq, r_init, p_init, m_init, curve, stepsize, condition_func)
    catch e
        @printf("error while solving diff equations: %s\n", e)
        #TODO: repeated code?
        curve.tvalues = curve.tvalues*LENGTH_UNIT_TO_SI*1e-3
        curve.xvalues = curve.xvalues
        curve.yvalues = curve.yvalues
        return curve
    end

    curve.tvalues = curve.tvalues*LENGTH_UNIT_TO_SI*1e-3
    curve.xvalues = curve.xvalues
    curve.yvalues = curve.yvalues

    return curve
end
