module TOV

#independent files
include("constants.jl")
include("diff.jl")
include("util.jl")

#dependent files
include("solver.jl")
include("solveutil.jl")

#Polytrope EoS modules
module WhiteDwarfPolytrope
    include("eos/polytropic_whitedwarf.jl")
end

module NeutronStarPolytrope
    include("eos/polytropic_neutronstar.jl")
end

end
