include("../constants.jl")

using .SIUnits

global const γ_rel = 4/3
global const K_REL = (ħ*c)/(12π^2)*((3π^2*Z)/(A*mₙ*c^2))^(γ_rel)

global const γ_nonrel = 5/3
global const K_NONREL = (ħ^2)/(15π^2*mₑ)*((3π^2*Z)/(A*mₙ*c^2))^(γ_nonrel)

function ϵ_rel(p::Real)::Real
    if p < 0 return 0 end
    return (p^(1/γ_rel))/(K_REL^(1/γ_rel))
end

function ϵ_nonrel(p::Real)::Real
    if p < 0 return 0 end
    return p^(1/γ_nonrel)
end
