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

# -------------------------------------------------------------------------
#
#   NEW CODE FOR SOLVING USING DIFFEQ LIBRARY
#
# -------------------------------------------------------------------------

using TOV
#NOTE: using DifferentialEquations package I should cite the author of this code
using SciMLBase
using OrdinaryDiffEq

const ErrorPressureZero = "error: pressure is lesser or equal to zero"

struct TOVSolution
    p₀::Real,
    r::AbstractVector{Real},
    p::AbstractVector{Real},
    M::AbstractVector{Real},
end

#NOTE r_max can be changed (for example in solving for white dwarfs)
function solve_tov(p₀::Real, eos::Function ; rinit::Real=0, rmax::Real=10e5*SI_TO_LENGTH_UNIT)::TOVSolution
    # don't know why I was still using 1e-24
    m_init = 0.0
    p_init = p₀

    pressure_eq(r, p, M) = begin
        r == 0 ? 0 : (-eos(p)*M/r^2)*(1 + p/eos(p))*(1 + 4πr^3*p/M)*(1 - 2M/r)^(-1)
    end
    mass_eq(r, p, M) = begin
        4π*r^2*eos(p)
    end

    f(du, u, p, t) = begin
        du[1] = pressure_eq(t, u[1], u[2])
        du[2] = mass_eq(t, u[1], u[2])
    end

    u0 = [p_init, m_init]
    tspan = (r_init, r_max)
    prob = ODEProblem(f, u0, tspan)

    # TODO: find a way to call a solve method that allows to break the solving at any point

    condition(u, t, integrator) = u[1] <= 0
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition, affect!)
    # Canonical Runge-Kutta Order 4 method. Uses adaptive stepping.
    sol = solve(prob, RK4(), callback = cb)
    return TOVSolution(p₀, sol.t, sol[1,:], sol[2,:])
end

function solve_tov(p₀::Real, eos::EOS ; rinit::Real=0, rmax::Real=10e5*SI_TO_LENGTH_UNIT)::TOVSolution
    return solve_tov(p₀,, eos.eos_fn, rinit=rinit, rmax=rmax)
end
