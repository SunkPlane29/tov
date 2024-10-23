#I will define some units in which c = G = M⊙ = 1, i will also define convertion factors for
#the basics units in SI: m, kg, s
#
#1(unit of mass)    ≈ 1.989*10^30kg
#1(unit of length)  ≈ 1.477*10^3m
#1(unit of time)    ≈ 4.927*10^(-6)s

function pressurediffeq(r::Real, P::Real, M::Real, ϵ::Function)::Real
    if r == 0
        return zero(r) # return 0 of type of r, this should make dynamic dispatch work better   
    end
    
    if P < 0
    	P = 0
    end

    -(ϵ(P)*M/r^2)*(1 + P/ϵ(P))*(1 + 4π*r^3*P/M)*(1 - 2M/r)^(-1)
end

function massdiffeq(r::Real, P::Real, M::Real, ϵ::Function)::Real
    if r == 0
        return zero(r)
    end
    
    if P < 0
    	P = 0
    end

    4π*r^2*ϵ(P)
end

function solvetov(P0::Real, ϵ::Function, h::Real=1m)::AbstractMatrix
    r0 = 1e-8
    m0 = 1e-24

    f(t, x) = [pressurediffeq(t, x[1], x[2], ϵ), massdiffeq(t, x[1], x[2], ϵ)]
    terminate(i, t, x) = x[1] <= 0 || x[1] < 1e-16

    sol = solvesystem(f, r0, [P0, m0], h, terminate)
    sol[:, 1] = sol[:, 1]*LENGTH_UNIT_TO_SI*1e-3
    sol[:, 2] = sol[:, 2]*PRESSURE_UNIT_TO_SI*JOULE_TO_MEV4*MEV4_TO_MEVFM3
    
    sol
end

function solvemrdiagram(P0::AbstractVector, ϵ::Function, h::Real=1m)::AbstractMatrix
    n = length(P0) 

    mrdiagram = zeros(n, 3)
    mrdiagram[:, 1] = P0
    mrdiagram[:, 2] = fill(0.0, n)
    mrdiagram[:, 3] = fill(0.0, n)

    Threads.@threads for i in 1:n
        sol = solvetov(P0[i], ϵ, h)
        mrdiagram[i, 2] = sol[end, 3]
        mrdiagram[i, 3] = sol[end, 1]
    end

    mrdiagram
end
