#I will define some units in which c = G = M⊙ = 1, i will also define convertion factors for
#the basics units in SI: m, kg, s
#
#1(unit of mass)    ≈ 1.989*10^30kg
#1(unit of length)  ≈ 1.477*10^3m
#1(unit of time)    ≈ 4.927*10^(-6)s

function pressure_diffeq(r::Real, p::Real, M::Real, ϵ::Function)::Real
    if r == 0
        return zero(r) # return 0 of type of r, this should make dynamic dispatch work better   
    end

    -(ϵ(p)*M/r^2)*(1 + p/ϵ(p))*(1 + 4π*r^3*p/M)*(1 - 2M/r)^(-1)
end

function mass_diffeq(r::Real, p::Real, M::Real, ϵ::Function)::Real
    if r == 0
        return zero(r)
    end

    4π*r^2*ϵ(p)
end

function solve_tov(p0::Real, ϵ::Function, h::Real=1m)::AbstractMatrix
    r0 = 1e-8
    m0 = 1e-24

    f(t, x) = [pressure_diffeq(t, x[1], x[2], ϵ), mass_diffeq(t, x[1], x[2], ϵ)]
    condition(i, t, x) = x[1] > 0

    sol = solve_system(f, r0, [p0, m0], h, condition)
    sol[:, 1] = sol[:, 1]*LENGTH_UNIT_TO_SI*1e-3
    sol[:, 2] = sol[:, 2]*PRESSURE_UNIT_TO_SI*JOULE_TO_MEV4*MEV4_TO_MEVFM3
    
    sol
end

function solve_mrdiagram(p0::AbstractVector, ϵ::Function, h::Real=1m)::AbstractMatrix
    n = length(p0) 

    mrdiagram = zeros(n, 3)
    mrdiagram[:, 1] = p0
    mrdiagram[:, 2] = fill(0.0, n)
    mrdiagram[:, 3] = fill(0.0, n)

    Threads.@threads for i in 1:n
        sol = solve_tov(p0[i], ϵ, h)
        mrdiagram[i, 2] = sol[end, 2]
        mrdiagram[i, 3] = sol[end, 1]
    end

    mrdiagram
end