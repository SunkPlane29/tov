using LinearAlgebra

struct LinearInterpolation
    x::Vector{Float64}
    y::Vector{Float64}
    a::Vector{Float64}
    b::Vector{Float64}
end

function LinearInterpolation(x::AbstractVector, y::AbstractVector)::LinearInterpolation
    n = length(x)-1
    a = zeros(n)
    b = zeros(n)

    for i in 1:n
        a[i] = (y[i+1] - y[i])/(x[i+1] - x[i])
        b[i] = y[i] - a[i]*x[i]
    end

    LinearInterpolation(x, y, a, b)
end

function (li::LinearInterpolation)(x::Real)::Real
    #extrapolation
    if x < li.x[1]
        return li.a[1]*x + li.b[1]
    elseif x > li.x[end]
        return li.a[end]*x + li.b[end]
    end

    i = searchsortedlast(li.x, x)

    if x == li.x[i]
        return li.y[i]
    end

    li.a[i]*x + li.b[i]
end

#TODO: implement derivative
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

struct EoS
    P::Vector{Float64}
    ϵ::Vector{Float64}
    itp::Union{LinearInterpolation, CubicSpline}
end

using CSV, DataFrames

#TODO: implement support for datafiles (not CSV) and other unit systems (not MeVfm3)
function EoS(file::AbstractString, header=["P", "ϵ"], method::Symbol=:cubic_spline)::EoS
    df = CSV.File(file, header=header) |> DataFrame
    if method == :linear_interpolation
        return EoS(df.P, df.ϵ, LinearInterpolation(df.P, df.ϵ))
    elseif method == :cubic_spline
        return EoS(df.P, df.ϵ, CubicSplineInterpolation(df.P, df.ϵ))
    end
end

function (eos::EoS)(P::Real)::Real
    eos.itp(P)
end