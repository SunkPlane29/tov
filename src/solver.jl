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
using FiniteDiff

struct TOVSolution
    p₀::Real
    r::AbstractVector
    p::AbstractVector
    m::AbstractVector
    phi::AbstractVector
    H::AbstractVector
    beta::AbstractVector
    k2::Real
    Lambda::Real
end

using Printf

function tov_eq(r::Real, p::Real, m::Real, phi::Real, H::Real, beta::Real, eos::Function)::Real
    -(eos(p)*m/r^2)*(1 + p/eos(p))*(1 + 4π*r^3*p/m)*(1 - 2m/r)^(-1)
end

function mass_continuity_eq(r::Real, p::Real, m::Real, phi::Real, H::Real, beta::Real, eos::Function)::Real
    4π*r^2*eos(p)
end

function phi_eq(r::Real, p::Real, m::Real, phi::Real, H::Real, beta::Real, eos::Function)::Real
    -tov_eq(r, p, m, phi, H, beta, eos)/(eos(p) + p)
end

function H_eq(r::Real, p::Real, m::Real, phi::Real, H::Real, beta::Real, eos::Function)::Real
    beta
end

#NOTE: might have an error in order of calculation
function beta_eq(r::Real, p::Real, m::Real, phi::Real, H::Real, beta::Real, eos::Function, eos_rev::Function)::Real
    f(p) = FiniteDiff.finite_difference_derivative(eos_rev, p)
    dbetadr = 2*(1 - 2m/r)^(-1)*H*(-2π*(5eos(p) + 9p + f(eos(p) + p)) + 3/r^2 +
        2*(1 - 2m/r)^(-1)*(m/r^2 + 4π*r*p)^2) +
        (2beta/r)*(1 - 2m/r)^(-1)*(-1 + m/r + 2π*r^2*(eos(p) - p))
    return dbetadr
end

#NOTE: for most stars, ϵ_sup=0
function k2(R::Real, M::Real, H::Real, beta::Real, ϵ_sup::Real)::Real
    y = R*beta/H - 4π*R^3*ϵ_sup/M
    C = M/R

    k_2 = (8C^(5)/5)*(1 - 2C)^(2)*(2 + 2C*(y - 1) - y)*(2C*(6 - 3y + 3C*(5y - 8)) +
        4C^(3)*(13 - 11y + C*(3y - 2) + 2C^(2)*(1 + y)) +
        3*(1 - 2C)^(2)*(2 - y + 2C*(y - 1))*log(1 - 2C))^(-1)

    return k_2
end

function Lambda(R::Real, M::Real, k2::Real)::Real
    C = M/R

    Λ = (2/3)*k2*C^(-5)

    return Λ
end

function get_diff_eq_system(eos::Function, eos_rev::Function)::Function
    dp(r, p, m, phi, H, beta) = begin
        tov_eq(r, p, m, phi, H, beta, eos)
    end
    dm(r, p, m, phi, H, beta) = begin
        mass_continuity_eq(r, p, m, phi, H, beta, eos)
    end
    dphi(r, p, m, phi, H, beta) = begin
        phi_eq(r, p, m, phi, H, beta, eos)
    end
    dH(r, p, m, phi, H, beta) = begin
        H_eq(r, p, m, phi, H, beta, eos)
    end
    dbeta(r, p, m, phi, H, beta) = begin
        beta_eq(r, p, m, phi, H, beta, eos, eos_rev)
    end

    f(du, u, p, t) = begin
        du[1] = dp(t, u[1], u[2], u[3], u[4], u[5])
        du[2] = dm(t, u[1], u[2], u[3], u[4], u[5])
        du[3] = dphi(t, u[1], u[2], u[3], u[4], u[5])
        du[4] = dH(t, u[1], u[2], u[3], u[4], u[5])
        du[5] = dbeta(t, u[1], u[2], u[3], u[4], u[5])
    end

    return f
end

function condition_func()
    condition(u, t, integrator) = begin
        if u[1] <= 0
            return true
        end

        return false
    end

    return condition
end

function solve_tov(p₀::Real, eos::Function, eos_rev::Function ; rmax::Real=10e5*SI_TO_LENGTH_UNIT,
                   stepsize::Real=1*SI_TO_LENGTH_UNIT, ϵsup::Real=0.0)::TOVSolution
    rinit = 1e-8
    pinit = p₀
    minit = 1e-24
    phiinit = 0
    Hinit = rinit^2
    betainit = 2*rinit

    f = get_diff_eq_system(eos, eos_rev)

    u0 = [pinit, minit, phiinit, Hinit, betainit]
    tspan = (rinit, rmax)
    prob = ODEProblem(f, u0, tspan)

    condition = condition_func()
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition, affect!)
    # Canonical Runge-Kutta Order 4 method. Uses adaptive stepping.
    sol = solve(prob, RK4(), callback = cb,
                adaptive=false, dt=stepsize)

    r = sol.t.*LENGTH_UNIT_TO_SI.*1e-3
    m = sol[2,:]
    H = sol[4,:]
    beta = sol[5,:]

    k_2 = k2(last(r), last(m), last(H), last(beta), ϵsup)
    Λ = Lambda(last(r), last(m), k_2)

    #NOTE: I have no idea of the units of ϕ, H and β so I will just output them in the units
    #of this system
    return TOVSolution(p₀, sol.t.*LENGTH_UNIT_TO_SI.*1e-3, sol[1,:]*PRESSURE_UNIT_TO_MEVFM3,
        sol[2,:], sol[3,:], sol[4,:], sol[5,:], k_2, Λ)
end

function solve_tov(p₀::Real, eos::EOS ; rmax::Real=10e5*SI_TO_LENGTH_UNIT,
                   stepsize::Real=1*SI_TO_LENGTH_UNIT, ϵsup::Real=0.0)::TOVSolution
    return solve_tov(p₀, eos.eos_function, eos.eos_function_rev, stepsize=stepsize, ϵsup=ϵsup)
end

struct SequenceSolution
    p₀::AbstractVector
    R::AbstractVector
    M::AbstractVector
    k2::AbstractVector
    Lambda::AbstractVector
    Phi::AbstractVector
    H::AbstractVector
    Beta::AbstractVector
end

function solve_sequence(p₀::AbstractVector, eos::Function, eos_rev::Function ;
                        rmax::Real=10e5*SI_TO_LENGTH_UNIT, stepsize::Real=1*SI_TO_LENGTH_UNIT,
                        ϵsup::Real=0.0)::SequenceSolution
    rinit = 1e-8
    pinit = p₀
    minit = 1e-24
    phiinit = 0
    Hinit = rinit^2
    betainit = 2*rinit

    f = get_diff_eq_system(eos, eos_rev)

    u0 = [pinit[1], minit, phiinit, Hinit, betainit]
    tspan = (rinit, rmax)
    prob = ODEProblem(f, u0, tspan)

    condition = condition_func()
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition, affect!)

    prob_func = (prob, i, repeat) -> remake(prob, u0=[pinit[i], minit, phiinit, Hinit, betainit])
    ensemble_prob = EnsembleProblem(prob, prob_func=prob_func, safetycopy=false)

    sol = solve(ensemble_prob, RK4(), EnsembleThreads(), callback=cb, trajectories=length(pinit),
                adaptive=false, dt=stepsize)

    R = []
    M = []
    k_2 = []
    Λ = []
    Phi = []
    H = []
    Beta = []

    for soli in sol
        Ri = last(soli.t)
        Mi = last(soli.u)[2]
        Hi = last(soli.u)[4]
        betai = last(soli.u)[5]
        phii = last(soli.u)[3]

        k_2i = k2(Ri, Mi, Hi, betai, ϵsup)
        Λi = Lambda(Ri, Mi, k_2i)

        append!(R, Ri*LENGTH_UNIT_TO_SI*1e-3)
        append!(M, Mi)
        append!(k_2, k_2i)
        append!(Λ, Λi)
        append!(Phi, phii)
        append!(H, Hi)
        append!(Beta, betai)
    end

    return SequenceSolution(p₀ .* PRESSURE_UNIT_TO_MEVFM3, R, M, k_2, Λ, Phi, H, Beta)
end

function solve_sequence(p₀::AbstractVector, eos::EOS ;
                        rmax::Real=10e5*SI_TO_LENGTH_UNIT, stepsize::Real=1*SI_TO_LENGTH_UNIT,
                        ϵsup::Real=0.0)::SequenceSolution
    return solve_sequence(p₀, eos.eos_function, eos.eos_function_rev, stepsize=stepsize, ϵsup=ϵsup)
end
