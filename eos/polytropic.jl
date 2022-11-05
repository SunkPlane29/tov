include("../constants.jl")

using .CGSUnits

global const γ_rel = 4/3
global const K_REL = (ħ*c)/(12π^2)*((3π^2*Z)/(A*mₙ*c^2))^(γ_rel)

global const γ_nonrel = 5/3
global const K_NONREL = (ħ^2)/(15π^2*mₑ)*((3π^2*Z)/(A*mₙ*c^2))^(γ_nonrel)

#TODO: define dimensionless ϵ using Lane-Emdem equation, search about this equation

function ϵ_rel(p::Real)::Real
    if p < 0 return 0 end
    return (p^(1/γ_rel))/(K_REL^(1/γ_rel))
end

function ϵ_nonrel(p::Real)::Real
    if p < 0 return 0 end
    return p^(1/γ_nonrel)
end

#TODO: plugging the same values in the same equation gives us different values from the tov_undergrad
function ϵ₀(α::Real)::Real
    R₀ = G*MSOLAR/c^2
    return ((1/K_REL)*(R₀/α)^γ_rel)^(1-γ_rel)
end
