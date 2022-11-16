include("../main.jl")
function main()
    let p₀ = 1.54e-16
        ϵ₀ = 5.61970127e+38
        r₀ = 0.9319e4

        if length(ARGS) == 4
            p₀ = parse(Float64, ARGS[2])
            ϵ₀ = parse(Float64, ARGS[3])
            r₀ = parse(Float64, ARGS[4])
        end

        solve_plot(p₀, ϵ₀, r₀)
    end
end
main()