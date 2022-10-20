include("eos.jl")
include("diff.jl")

SOLAR_MASS = 1.241105355e43      # MeV
Λ = 587.9                        # MeV
G = 2.44/(Λ^2)                   # in MeV^-2
c = 1                            # dimensionless
α = -(G * SOLAR_MASS)/c^2        # kg (i suppose)
β = (4π)/(SOLAR_MASS * c^2)      # kg^-1 (i suppose)
m₀ = 5.6                         # MeV

# -------------------------------------------------------
#
# PRIMEIRA TENTATIVA (acho que não é a primeira de verdade)
#
#   - agora que eu parei pra pensar, talvez eu tenha que pegar a equação da reta da EOS
#   pra ter uma valor pra todo p. Se eu ficar no arquivo de dados, talvez não dê tão certo
#   na hora de resolver as equações acopladas (não vai casar bem com os steps que são feitos
#   no método de resolução)
#
# -------------------------------------------------------

eos = loadEOS()
avg_slope = average_slope(eos)
ϵ₀ = eos.df.energy_density[1]

# convert from fm^-4 to MeV^4
# never utilized this hehe
function pfm_to_pmev(fm::Real)::Real
    return fm * 1.5161165184e9
end

# pressure differential equation with some change for using
function pressure_equation(r::Real, p::Real, M::Real)::Real
    # TODO: make function to integrate M(r) and one functoin to get values of energy density
    # corresponding to given pressure from the data file

    if r == 0
        throw("r = 0 means division by 0 in pressure equation")
    end

    if M == 0
        throw("M = 0 means a 0 angular coefficient for the tangent, that can possibly cause the solving to stay at 0 throughout")
    end

    return α*linearget_energydensity(avg_slope, ϵ₀, p)*mass_function(r, p, M)/r^2
end

# is the actual coupled differential equation
function mass_equation(r::Real, p::Real, M::Real)::Real
    return β * r^2*linearget_energydensity(avg_slope, ϵ₀, p)
end

using QuadGK
# returns a mass value corresponding to the mass in the interior of the radius r
function mass_function(r::Real, p::Real, M::Real)::Real
    mm(r) = r^2*linearget_energydensity(avg_slope, ϵ₀, p)

    (int, err) = quadgk(mm, 0, r)
    return β * int
end

function solve_tov()::Tuple{Curve,Curve}
    p₀ = last(eos.df.pressure)
    r₀ = 1e-3
    stepsize = 0.08
    n = 100

    return solve_system(pressure_equation, mass_equation, r₀, p₀, m₀, stepsize, n)
end
