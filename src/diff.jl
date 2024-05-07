function next_point(f::Function, t::Real, x::AbstractVector, h::Real)::AbstractVector
    kn1 = f(t, x)
    kn2 = f(t + h/2, x .+ (h/2).*kn1)
    kn3 = f(t + h/2, x .+ (h/2).*kn2)
    kn4 = f(t + h, x .+ h.*kn3)

    x .+ (h/6).*(kn1 .+ 2kn2 .+ 2kn3 .+ kn4)
end

function solve_system(f::Function, t0::Real, x0::AbstractVector, h::Real, condition::Function, maxiter::Integer=100000)::AbstractMatrix
    tn = t0
    xn = x0
    sol = [[t0, x0...]]
    sizehint!(sol, maxiter)
    i = 1

    while condition(i, tn, xn) && i <= maxiter
        xn = next_point(f, tn, xn, h)
        tn += h
        push!(sol, [tn, xn...])
        i += 1
    end

    permutedims(hcat(sol...))
end