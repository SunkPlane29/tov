include("../constants.jl") 
using .NaturalUnits

global const γ_rel = 4/3
global const K_REL = (ħ*c)/(12π^2)*((3π^2*Z)/(A*mₙ*c^2))^(γ_rel)

global const γ_nonrel = 5/3
global const K_NONREL = (ħ^2)/(15π^2*mₑ)*((3π^2*Z)/(A*mₙ*c^2))^(γ_nonrel)

#TODO: do equation of state with dimensionless p and ϵ defined using a central energy density
#while also defining dimensionless M that corresponds to solar masses and finally define
#dimensionless r that corresponds to a fraction of the total radius of the star (this one I'm not sure)
#the plan is defining all the using dimensionless values, and, knowing the convertions, convert back all
#the results in the desired dimensions

function ϵ_rel(p::Real)::Real
    if p < 0 return 0 end
    return (p^(1/γ_rel))/(K_REL^(1/γ_rel))
end

function ϵ_nonrel(p::Real)::Real
    if p < 0 return 0 end
    return p^(1/γ_nonrel)
end

using Printf

#TODO: plugging the same values in the same equation gives us different values from the tov_undergrad
function ϵ₀(α::Real)::Real
    R₀ = G*MSOLAR/c^2
    e0 = ((1/K_REL)*(R₀/α)^γ_rel)^(1-γ_rel)
    β = (4π*e0)/(MSOLAR*c^2*(K_REL*e0^(1-γ_rel))^(1/γ_rel))
    @printf("with α = %.8e, fixed ϵ₀ = %.8e and β = %.8e\n", α, e0, β)
    return e0
end
