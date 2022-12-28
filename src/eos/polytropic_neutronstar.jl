global const A::Real = 56
global const Z::Real = 26

global const γ_nonrel::Real = 5/3
global const K_NONREL::Real = (ħ^2)/(15π^2*mₙ)*((3π^2*Z)/(A*mₙ*c^2))^(γ_nonrel)
global const γ_rel::Real = 1
global const K_REL::Real = 1/3

#polytrope equation of state that returns the energy density given a pressure. γ and K
#are just constants that differentiate the relativistic and non-relativistic limits
function polytrope(p::Real, γ::Real, K::Real)::Real
    return (p^(1/γ))/(K^(1/γ))
end
