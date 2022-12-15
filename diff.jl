mutable struct Curve
    tvalues::Vector{Float64}
    xvalues::Vector{Float64}
    yvalues::Vector{Float64}
end

#TODO: there is a way to do this as a vector
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
