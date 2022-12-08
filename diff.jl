mutable struct Curve
    tvalues::AbstractVector{Real}
    xvalues::AbstractVector{Real}
    yvalues::AbstractVector{Real}
end

#TODO: there is a way to do this as a vector
function next_point(f::Function, g::Function, t::Real, x::Real, y::Real, stepsize::Real)::Tuple{Real,Real,Real}
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

    return (next_t, next_x, next_y)
end

#FIXME: there is one big big big problem in this method, while it does match the curve pretty well, changes in the
#stepsize are affecting the solution, this should not be the case as we want the stepsize to get smaller and the
#curves more precise. Find the cause of this issue
function solve_system!(x::Function, y::Function, t₀::Real, x₀::Real, y₀::Real, curve::Curve, stepsize::Real, condition::Function)
    append!(curve.tvalues, t₀)
    append!(curve.xvalues, x₀)
    append!(curve.yvalues, y₀)

    previoust = t₀
    previousx = x₀
    previousy = y₀
    nt, nx, ny = 1, 1, 1

    i = 1

    while condition(i, previoust, previousx, previousy)
        (nt, nx, ny) = next_point(x, y, previoust, previousx, previousy, stepsize)

        append!(curve.tvalues, nt)
        append!(curve.xvalues, nx)
        append!(curve.yvalues, ny)

        previoust = nt
        previousx = nx
        previousy = ny
        i += 1
    end
end
