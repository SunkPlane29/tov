include("main.jl")

function ajust_parameters(p₀::Real, desired_radius::Real, desired_mass::Real, γ::Real, K::Real, radius_precision::Real = 1000, mass_precision::Real = 0.05)
    ϵ₀_min = 1.0e32
    ϵ₀_max = 1.0e40
    r₀_min = 1e3
    r₀_max = 1e5
    n = 1000
    ϵ₀_del = (ϵ₀_max-ϵ₀_min)/n
    r₀_del = (r₀_max-r₀_min)/n

    ϵ₀_min_diff = 0
    r₀_min_diff = 0
    R_min_diff = 1e10
    M_min_diff = 1e10

    let ϵ₀, r₀
        for i in 1:n
            r₀ = r₀_min+(i-1)*r₀_del
            for j in 1:n
                ϵ₀ = ϵ₀_min+(j-1)*ϵ₀_del

                curve = solve(p₀, ϵ₀, r₀, γ, K, false)

                R_diff = abs(last(curve.tvalues) - desired_radius)
                M_diff = abs(last(curve.yvalues) - desired_mass)
                if R_diff < R_min_diff && M_diff < M_min_diff
                    R_min_diff = R_diff
                    ϵ₀_min_diff = ϵ₀
                    r₀_min_diff = r₀
                end


                if R_diff < radius_precision && M_diff < mass_precision
                    @printf("found R = %.8e and M = %.8e with parameters ϵ₀ = %.8e and r₀ = %.8e\n", last(curve.tvalues), last(curve.yvalues), ϵ₀, r₀)
                    return
                end
            end
        end
    end

    @printf("no R and M in desired accuracy, found only ϵ₀ = %.8e and r₀ = %.8e with a difference of R %.8e and difference of M %.8e\n", ϵ₀_min_diff, r₀_min_diff, R_min_diff, M_min_diff)
end
