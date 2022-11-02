struct Curve
    tvalues::AbstractVector{Real}
    xvalues::AbstractVector{Real}
    yvalues::AbstractVector{Real}
end

#TODO: there is a way to do this as a vector
function next_point(f::Function, g::Function, t::Real, x::Real, y::Real, stepsize::Real)::Tuple{Real,Real,Real}
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

    return (next_t, next_x, next_y)
end

#TODO: maybe make this accept a function (with some arguments, like x, y, n) that determines if the loop should stop
#maybe make this function in another julia method (remember, a julia method is a function)
function solve_system(x::Function, y::Function, t₀::Real, x₀::Real, y₀::Real, stepsize::Real, n::Integer)::Curve
    curve = Curve(Real[t₀], Real[x₀], Real[y₀])

    previoust = t₀
    previousx = x₀
    previousy = y₀
    nt, nx, ny = 1, 1, 1

    i = 1
    eps = 1.0e-17

    while nx > 0
    #while i <= n
    #while abs(nx) > eps
        #println(i)
        #if abs(nx) < eps break end
        if i > n break end
        (nt, nx, ny) = next_point(x, y, previoust, previousx, previousy, stepsize)

        append!(curve.tvalues, nt)
        append!(curve.xvalues, nx)
        append!(curve.yvalues, ny)

        previoust = nt
        previousx = nx
        previousy = ny
        i += 1
    end

    return curve
end
