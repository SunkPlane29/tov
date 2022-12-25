module TOV

#independent files
include("constants.jl")
export c, G, ħ, mₙ, mₑ, si_c, si_G, si_MSOLAR, si_ħ, si_mₙ, si_mₑ, MSOLAR, MASS_UNIT_TO_SI, LENGTH_UNIT_TO_SI,
    TIME_UNIT_TO_SI, SI_TO_MASS_UNIT, SI_TO_LENGTH_UNIT, SI_TO_TIME_UNIT, PRESSURE_UNIT_TO_SI,
    SI_TO_PRESSURE_UNIT, SI_TO_GEV_FM3, GEV_FM3_TO_SI, ħc, PRESSURE_FM_TO_MEV, SI_EV_TO_JOULE, MEV4_TO_JOULE

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
    include("constants.jl")
    include("eos/polytropic_whitedwarf.jl")
    export A, Z, γ_nonrel, K_NONREL, γ_rel, K_REL, polytrope
end

module NeutronStarPolytrope
    include("constants.jl")
    include("eos/polytropic_neutronstar.jl")
    export A, Z, γ_nonrel, K_NONREL, γ_rel, K_REL, polytrope
end

end
