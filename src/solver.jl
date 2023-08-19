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

# using TOV
#NOTE: using DifferentialEquations package I should cite the author of this code
using SciMLBase
using OrdinaryDiffEq

const ErrorPressureZero = "error: pressure is lesser or equal to zero"

struct TOVSolution
    p₀::Real
    r::AbstractVector{Real}
    p::AbstractVector{Real}
    M::AbstractVector{Real}
end

using Printf

function tov_eq(r::Real, p::Real, M::Real, eos::Function)::Real
    -(eos(p)*M/r^2)*(1 + p/eos(p))*(1 + 4π*r^3*p/M)*(1 - 2M/r)^(-1)
end

function mass_continuity_eq(r::Real, p::Real, M::Real, eos::Function)::Real
    4π*r^2*eos(p)
end

function get_diff_eq_system(eos::Function)::Function
    pressure_eq(r, p, M) = begin
        tov_eq(r, p, M, eos)
    end
    mass_eq(r, p, M) = begin
        mass_continuity_eq(r, p, M, eos)
    end

    f(du, u, p, t) = begin
        du[1] = pressure_eq(t, u[1], u[2])
        du[2] = mass_eq(t, u[1], u[2])
    end

    return f
end

function solve_tov(p₀::Real, eos::Function ; rinit::Real=1e-8, minit::Real=1e-24, rmax::Real=10e5*SI_TO_LENGTH_UNIT, eps::Real=1e-10)::TOVSolution
    pinit = p₀
    if rinit == 0
        rinit = 1e-8
    end

    f = get_diff_eq_system(eos)

    u0 = [pinit, minit]
    tspan = (rinit, rmax)
    prob = ODEProblem(f, u0, tspan)

    rprev = 0.0
    mprev = 0.0
    pprev = 0.0
    #first 10 or so could trick these into thinking its converging
    i = 0
    condition(u, t, integrator) = begin
        rtemp = rprev
        ptemp = pprev
        mtemp = mprev
        rprev = t
        pprev = u[1]
        mprev = u[2]
        i += 1

        #NOTE: pressure convergence appears to be better
        if rtemp == 0.0
            return false
        elseif abs(ptemp - u[1]) <= eps && i > 10
            return true
        end

        return false
    end
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition, affect!)
    # Canonical Runge-Kutta Order 4 method. Uses adaptive stepping.
    sol = solve(prob, RK4(), callback = cb)
    return TOVSolution(p₀, sol.t.*LENGTH_UNIT_TO_SI.*1e-3, sol[1,:]*PRESSURE_UNIT_TO_MEVFM3, sol[2,:])
end

function solve_tov(p₀::Real, eos::EOS ; rinit::Real=1e-8, minit::Real=1e-24, rmax::Real=10e5*SI_TO_LENGTH_UNIT, eps::Real=1e-10)::TOVSolution
    return solve_tov(p₀, eos.eos_function, rinit=rinit, rmax=rmax, eps=eps)
end

struct SequenceSolution
    p₀::AbstractVector{Real}
    R::AbstractVector{Real}
    M::AbstractVector{Real}
end

function solve_sequence(p₀::AbstractVector{Real}, eos::Function ; rinit::Real=0.0, minit::Real=1e-24, rmax::Real=10e5*SI_TO_LENGTH_UNIT, eps::Real=1e-10)
    minit = fill(1e-24, length(p₀))
    pinit = p₀

    f = get_diff_eq_system(eos)

    u0 = [pinit[1], minit[1]]
    tspan = (rinit, rmax)
    prob = ODEProblem(f, u0, tspan)

    #FIXME: ugly repeated code here, don't know how to make it better
    rprev = 0.0
    mprev = 0.0
    pprev = 0.0
    i = 0
    condition(u, t, integrator) = begin
        rtemp = rprev
        ptemp = pprev
        mtemp = mprev
        rprev = t
        pprev = u[1]
        mprev = u[2]
        i += 1

        if rtemp == 0.0
            return false
        elseif abs(ptemp - u[1]) <= eps && i > 10
            return true
        end

        return false
    end
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition, affect!)

    prob_func = (prob, i, repeat) -> remake(prob, u0=[pinit[i], minit[1]])
    ensemble_prob = EnsembleProblem(prob, prob_func=prob_func, safetycopy=false)

    sol = solve(ensemble_prob, RK4(), EnsembleThreads(), trajectories=length(pinit))

    return sol
end

function solve_sequence(p₀::AbstractVector{Real}, eos::EOS ; rinit::Real=1e-8, minit::Real=1e-24, rmax::Real=10e5*SI_TO_LENGTH_UNIT, eps::Real=1e-10)
    return solve_sequence(p₀, eos.eos_function, rinit=rinit, minit=minit, rmax=rmax, eps=eps)
end
