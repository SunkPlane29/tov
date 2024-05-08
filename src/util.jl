#TODO: implement cubic spline
#TODO: implement derivative (of cubic spline if necessary)

using LinearAlgebra

struct CubicSpline
    x::Vector{Float64}
    y::Vector{Float64}
    k::Vector{Float64}
    c::Vector{Float64}
    d::Vector{Float64}
end

function CubicSplineInterpolation(x::AbstractVector, y::AbstractVector)::CubicSpline
    n = length(x)-1
    M = zeros(n+1, n+1)
    b = zeros(n+1)

    M[1, :] = [2/(x[2] - x[1]) ; 1/(x[2] - x[1]) ; zeros(n-1)]
    b[1] = 3*(y[2] - y[1])/(x[2] - x[1])^2

    for i in 2:n
        M[i, :] = [zeros(i-2) ; 1/(x[i] - x[i-1]) ; 2*(1/(x[i] - x[i-1]) + 1/(x[i+1] - x[i])) ; 1/(x[i+1] - x[i]) ; zeros(n-i)]
        b[i] = 3*((y[i] - y[i-1])/(x[i] - x[i-1])^2 + (y[i+1] - y[i])/(x[i+1] - x[i])^2)
    end

    M[n+1, :] = [zeros(n-1) ; 1/(x[n+1] - x[n]) ; 2/(x[n+1] - x[n])]
    b[n+1] = 3*(y[n+1] - y[n])/(x[n+1] - x[n])^2

    k = M\b
    c = zeros(n)
    d = zeros(n)
    for i in 1:n
        c[i] = k[i]*(x[i+1] - x[i]) - (y[i+1] - y[i])
        d[i] = -k[i+1]*(x[i+1] - x[i]) + (y[i+1] - y[i])
    end
    
    CubicSpline(x, y, k, c, d)
end

function (cs::CubicSpline)(x::Real)::Real
    if x < cs.x[1] || x > cs.x[end]
        throw(ArgumentError("x must be within the range of x"))
    end

    i = searchsortedlast(cs.x, x)

    if x == cs.x[i]
        return cs.y[i]
    end

    t = (x - cs.x[i])/(cs.x[i+1] - cs.x[i])

    (1 - t)*cs.y[i] + t*cs.y[i+1] + t*(1 - t)*(cs.c[i]*(1 - t) + cs.d[i]*t) 
end