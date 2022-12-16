include("../main.jl")
function main()
    let p₀ = 1.603e22

        if length(ARGS) == 4
            p₀ = parse(Float64, ARGS[2])
        end

        solve_plot(p₀, ϵ₀, r₀)
    end
end
main()
