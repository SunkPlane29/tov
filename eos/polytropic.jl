const A = 56
const Z = 26

global const γ_nonrel = 5/3
global const K_NONREL = (ħ^2)/(15π^2*mₙ)*((3π^2*Z)/(A*mₙ*c^2))^(γ_nonrel)

function polytrope(p::Real, γ::Real, K::Real)::Real
    return (p^(1/γ))/(K^(1/γ))
end
