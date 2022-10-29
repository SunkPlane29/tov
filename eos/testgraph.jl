using Plots

global const mₑ = 9.1093837e-28            # electron mass in gm
global const mₙ = 1.674927471e-24          # nucleon mass in gm
global const c = 29979245800               # cm/s
global const ħ = 1.05457182e-27            # ergs s
global const ϵ₀ = 1.440559724e24 # constant ϵ₀ in ergs/cm^3
global const AZRATIO = 2.15                # ratio of nucleons per electrons

function electron_numberdensity(k_F::Real)::Real
    return (k_F^3)/(3π^2 * ħ^3)
end

function electron_energydensity(k_F::Real)::Real
    x = k_F/(mₑ*c)
    return (ϵ₀/8)*((2x^3 + x)*sqrt(1+x^2) - asinh(x))
end

function energy_density(k_F::Real)::Real
    return (electron_numberdensity(k_F)*mₙ * AZRATIO) + electron_energydensity(k_F)
end

# ta dando algum erro aqui, não se no que
function pressure(k_F::Real)::Real
    x = k_F/(mₑ*c)
    return (ϵ₀/24)*((2x^3 - 3x)*sqrt(1+x^2) + 3asinh(x))
end

function main()
    evalues = []
    pvalues = []
    kfvalues = []

    kfame = 0
    kfbme = 2
    n = 1000
    stepsize = (kfbme - kfame)/(n)
    for i = 1:n
        kfvalues = append!(kfvalues, kfame + stepsize*i*mₑ)
    end

    for k_F in kfvalues
        append!(evalues, energy_density(k_F))
        append!(pvalues, pressure(k_F))
    end

    p = plot(evalues, pvalues)
    savefig(p, "plot.png")
end
