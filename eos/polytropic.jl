global const γ_rel = 4/3
global const K_REL = (ħ*c)/(12π^2)*((3π^2*Z)/(A*mₙ*c^2))^(γ_rel)

global const γ_nonrel = 5/3
global const K_NONREL = (ħ^2)/(15π^2*mₑ)*((3π^2*Z)/(A*mₙ*c^2))^(γ_nonrel)

function polytrope(p::Real, γ::Real, K::Real)::Real
    return (p^(1/γ))/(K^(1/γ))
end

function ϵ₀_const(ϵ₀::Real, γ::Real)::Real
    return 1/ϵ₀^((γ-1)/γ)
end
