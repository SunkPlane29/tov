function nextpoint(f::Function, t::Real, x::AbstractVector, h::Real)::AbstractVector
    kn1 = f(t, x)
    kn2 = f(t + h/2, x .+ (h/2).*kn1)
    kn3 = f(t + h/2, x .+ (h/2).*kn2)
    kn4 = f(t + h, x .+ h.*kn3)

    x .+ (h/6).*(kn1 .+ 2kn2 .+ 2kn3 .+ kn4)
end

function solvesystem(f::Function, t0::Real, x0::AbstractVector, h::Real, terminate::Function, maxiter::Integer=100000)::AbstractMatrix
    tn = t0
    xn = x0
    sol = [[t0, x0...]]
    sizehint!(sol, maxiter)
    i = 1

    for i in 1:maxiter 
        xn = nextpoint(f, tn, xn, h)
        tn += h
        if terminate(i, tn, xn)
            break
        end
        
        push!(sol, [tn, xn...])
    end

    permutedims(hcat(sol...))
end