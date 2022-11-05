#FIXME: no standard in decimals
#FIXME: these constants are ALL messed up, should use all of them in cgs standard and convert
#everything elso to this standard (i think i need to get everything to the same units)

global const Z      = 6                 #number of protons/electrons in a carbom atom
global const A      = 6                #number of neutrons in a carbon atom

module NaturalUnits
export G, c, ħ, mₑ, mₙ, MSOLAR

global const c      = 1                 #dimensionless
global const ħ      = 1                 #dimensionless
global const G      = 1.0e-51           #MeV⁻² (in natural units)
global const mₑ     = 0.51099895000     #electron mass in MeV (in natural units)
global const mₙ     = 939.56542052      #neutron mass in MeV (in natural units)
global const MSOLAR = 1.115449865e60    #solar mass in MeV
end

module SIUnits
export G, c, ħ, mₑ, mₙ, MSOLAR

global const G      = 6.674e-11         #m³⋅kg⁻¹⋅s⁻²
global const c      = 2.99792458e8      #m⋅s⁻²
global const ħ      = 1.0545919e-34     #J⋅s
global const mₑ     = 9.10938370e-31    #kg
global const mₙ     = 1.67492750e-27    #kg
global const MSOLAR = 1.98847e30        #kg
end

module CGSUnits
export G, c, ħ, mₑ, mₙ, MSOLAR

global const G      = 6.674e-8          #dyne⋅cm²⋅g⁻²
global const c      = 2.99792458e10     #cm⋅s⁻²
global const ħ      = 1.0545919e-27     #erg⋅s
global const mₑ     = 9.10938370e-28    #g
global const mₙ     = 1.67492750e-24    #g
global const MSOLAR = 1.98847e33        #g
end

module GeometrizedSolarUnits
export G, c, ħ, mₑ, mₙ, MSOLAR

global const G      = 1                                     #dimensionless
global const c      = 1                                     #dimensionless
global const MSOLAR = 1                                     #dimensionless
# I see some numerical instability comming from here, maybe this is wrong
global const ħ      = 1.0545919e-27 * 5.5953e-55 * 2.0296e5 #dimensionless
global const mₑ     = 9.10938370e-28 * 5.0279*10^-34        #dimensionless
global const mₙ     = 1.67492750e-24 * 5.0279*10^-34        #dimensionless
end

module GeometrizedUnits
export G, c, ħ, mₑ, mₙ, MSOLAR

global const G      = 1                                                 #dimensionless
global const c      = 1                                                 #dimensionless
# Probably another possible cause of numerical instability
global const ħ      = 1.0545919e-27 * 8.2627e-50 * 3.33564095198e-11    #cm²
global const mₑ     = 9.10938370e-28 * 7.4261e-29                       #cm
global const mₙ     = 1.67492750e-24 * 7.4261e-29                       #cm
global const MSOLAR = 1.98847e33 * 7.4261e-29                           #cm
end

# ------------------------------------------------------------
#
# CONVERSION FACTORS
#
# ------------------------------------------------------------

global const KM_MEV_convertionfac   = 6.241509074e12   #MeV⁻¹
global const MEV_KM_convertionfac   = 1.602176634e-13  #km
global const FM_MEV_4_convertionfac = 1.515333907e30   #MeV⁴
global const MEV_FM_4_convertionfac = 6.599205597e-31  #fm⁻⁴

global const bruno_p0 = 6.55264111653593790e-3*FM_MEV_4_convertionfac
