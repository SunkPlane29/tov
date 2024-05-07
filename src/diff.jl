#TODO: refactor

#Curve holds the solution data of two couppled differential equations, x and y,  that depend on a single
#variable t
mutable struct Curve
    tvalues::Vector{Float64}
    xvalues::Vector{Float64}
    yvalues::Vector{Float64}
end

#next_point calculates the next_point in the solution of the couppled differential equation x and y using Range-Kutta 4
#method
function next_point(f::Function, g::Function, t::Real, x::Real, y::Real, stepsize::Real)::Tuple{Real,Real,Real}
    h = stepsize
    v = [x,y]

    f_mod(t, v) = [f(t, v[1], v[2]), g(t, v[1], v[2])]

    kn1 = f_mod(t, v)
    kn2 = f_mod(t + h/2, v .+ (h/2)*kn1)
    kn3 = f_mod(t + h/2, v .+ (h/2)*kn2)
    kn4 = f_mod(t + h, v .+ h*kn3)

    next_v = v + ((h/6).*(kn1 .+ 2kn2 .+ 2kn3 .+ kn4))

    next_t = t + h
    next_x = next_v[1]
    next_y = next_v[2]

    return (next_t, next_x, next_y)
end

#solve_system! solves the system of two couppled differential equations x and y, that depend on a variable t. The solution
#of this system is stored in a Curve structure that should be passed as an argument. A condition function should also be
#passed as argument, this function will accept i (the index of the solution), previoust (previous value of t variable),
#previousx (previous value of the x variable) and prevousy (previous value of y variable)
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
