module TOV

#independent files
include("constants.jl")
export c, G, ħ, mₙ, mₑ, MSOLAR, MASS_UNIT_TO_SI, LENGTH_UNIT_TO_SI, TIME_UNIT_TO_SI, SI_TO_MASS_UNIT,
    SI_TO_LENGTH_UNIT, SI_TO_TIME_UNIT, PRESSURE_UNIT_TO_SI, SI_TO_PRESSURE_UNIT, SI_TO_GEV_FM3

include("diff.jl")
export Curve, next_point, solve_system!

include("util.jl")
export plot_curves, write_data, plot_from_datafile

#dependent files
include("solver.jl")
export solve_tov

include("solveutil.jl")
export solve, solve_plot, solve_data, solve_star_curve

#Polytrope EoS modules
module WhiteDwarfPolytrope
    #TODO: find where to put this
    function get_polytrope(γ::Real, K::Real)::Function
        eos(p) = polytrope(p, γ, K)
        return eos
    end

    include("constants.jl")
    include("eos/polytropic_whitedwarf.jl")
    export A, Z, γ_nonrel, K_NONREL, γ_rel, K_REL, polytrope, get_polytrope
end

module NeutronStarPolytrope
    function get_polytrope(γ::Real, K::Real)::Function
        eos(p) = polytrope(p, γ, K)
        return eos
    end

    include("constants.jl")
    include("eos/polytropic_neutronstar.jl")
    export A, Z, γ_nonrel, K_NONREL, γ_rel, K_REL, polytrope, get_polytrope
end

end
