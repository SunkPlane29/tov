include("util.jl")

struct Curve
    tvalues::AbstractVector{Real}
    yvalues::AbstractVector{Real}
end

struct Point
    t::Real
    x::Real
    y::Real
end

function next_point(f::Function, g::Function, t::Real, x::Real, y::Real, stepsize::Real)::Point
    #TODO: maybe use a better step size later

    #h = 1e-3
    h = stepsize

    kn1 = f(t, x, y)
    ln1 = g(t, x, y)

    kn2 = f(t + (1/2)*h, x + (1/2)*kn1, y + (1/2)*ln1)
    ln2 = g(t + (1/2)*h, x + (1/2)*kn1, y + (1/2)*ln1)

    kn3 = f(t + (1/2)*h, x + (1/2)*kn2, y + (1/2)*ln2)
    ln3 = g(t + (1/2)*h, x + (1/2)*kn2, y + (1/2)*ln2)

    kn4 = f(t + h, x + kn3, y + ln3)
    ln4 = g(t + h, x + kn3, y + ln3)

    kn = (1/6)*(kn1 + 2kn2 + 2kn3 + kn4)
    ln = (1/6)*(ln1 + 2ln2 + 2ln3 + ln4)

    next_t = t + h
    next_x = x + kn
    next_y = y + ln

    return Point(next_t, next_x, next_y)
end

#TODO: maybe make this accept a function (with some arguments, like x, y, n) that determines if the loop should stop
#maybe make this function in another julia method (remember, a julia method is a function)
function solve_system(x::Function, y::Function, t₀::Real, x₀::Real, y₀::Real, stepsize::Real, n::Integer)::Tuple{Curve,Curve}
    xcurve = Curve(Real[0], Real[x₀])
    ycurve = Curve(Real[0], Real[y₀])

    previoust = t₀
    previousx = x₀
    previousy = y₀

    i = 1

    p = Point(1, 1, 1)

    #TODO: change for a while loop (to integrate this while p ≠ 0)
    while p.x > 0
        #FIXME: this code is ugly
        if get_debug() println(i) end

        p = next_point(x, y, previoust, previousx, previousy, stepsize)
        t = p.t #any of these two works

        append!(xcurve.tvalues, t)
        append!(xcurve.yvalues, p.x)
        append!(ycurve.tvalues, t)
        append!(ycurve.yvalues, p.y)

        previoust = t
        previousx = p.x
        previousy = p.y
        i += 1
    end

    return (xcurve, ycurve)
end
