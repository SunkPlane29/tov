#mostly taken from CODATA (NIST)

global const c = G = MSOLAR = 1

global const si_c::Real         = 2.99792458e8      # m s⁻¹
global const si_G::Real         = 6.67430e-11       # m³ kg⁻¹ s⁻²
global const si_MSOLAR::Real    = 1.98855e30        # kg
global const si_ħ::Real         = 1.054571817e-34   # J s
global const si_mₙ::Real        = 1.67492749804e-27 # kg
global const si_mₑ::Real        = 9.1093837015e-31  #gk

global const MASS_UNIT_TO_SI::Real       = si_MSOLAR                # kg
global const LENGTH_UNIT_TO_SI::Real     = si_G*si_MSOLAR/si_c^2    # m
global const TIME_UNIT_TO_SI::Real       = si_G*si_MSOLAR/si_c^3    # s
global const SI_TO_MASS_UNIT::Real       = MASS_UNIT_TO_SI^(-1)     # dimensionless
global const SI_TO_LENGTH_UNIT::Real     = LENGTH_UNIT_TO_SI^(-1)   # dimensionless
global const SI_TO_TIME_UNIT::Real       = TIME_UNIT_TO_SI^(-1)     # dimensionless

#TODO: there may be issues with precision regarding this convertion factor
#also, I don't know how to handle uncertainties and errors, also decimals to use
global const PRESSURE_UNIT_TO_SI::Real = MASS_UNIT_TO_SI*LENGTH_UNIT_TO_SI^(-1)*TIME_UNIT_TO_SI^(-2) #Pa
global const SI_TO_PRESSURE_UNIT::Real = PRESSURE_UNIT_TO_SI^(-1)                                    #dimensionless

#don't know how I defined this but here it is
global const SI_TO_GEV_FM3::Real         = 6.241509074e−36    #GeV fm⁻³
global const GEV_GM3_TO_SI::Real         = SI_TO_GEV_FM3^(-1) #J m⁻³

global const ħ::Real  = si_ħ*SI_TO_MASS_UNIT*SI_TO_LENGTH_UNIT^2*SI_TO_TIME_UNIT^(-1) #dimensionless
global const mₙ::Real = si_mₙ*SI_TO_MASS_UNIT                                         #dimensionless
global const mₑ::Real = si_mₑ*SI_TO_MASS_UNIT                                         #dimensionless
