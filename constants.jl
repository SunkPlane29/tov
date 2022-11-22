#FIXME: no standard in decimals

global const Z      = 6                 #number of protons/electrons in a carbom atom
global const A      = 12                #number of mass in a carbon atom

module NaturalUnits
export G, c, ħ, mₑ, mₙ, MSOLAR, MEV_TO_KG_CONVERTIONFAC, LENGTH_TO_KM_CONVERTIONFACTOR

global const MEV_TO_KG_CONVERTIONFAC = 1.0e30/1.782661922
global const LENGTH_TO_KM_CONVERTIONFACTOR = 1.9733e-16          #MeV to km

global const c      = 1                                     #dimensionless
global const ħ      = 1                                     #dimensionless
global const G      = 6.67259e-45*197.327                   #fm⋅MeV⁻¹ (in natural units)
global const mₑ     = 0.51099895000                         #electron mass in MeV (in natural units)
global const mₙ     = 939.56542052                          #neutron mass in MeV (in natural units)
global const MSOLAR = 1.98847e30*MEV_TO_KG_CONVERTIONFAC    #solar mass in MeV
end

module SIUnits
export G, c, ħ, mₑ, mₙ, MSOLAR, LENGTH_TO_KM_CONVERTIONFACTOR

global const G      = 6.674e-11         #m³⋅kg⁻¹⋅s⁻²
global const c      = 2.99792458e8      #m⋅s⁻²
global const ħ      = 1.0545919e-34     #J⋅s
global const mₑ     = 9.10938370e-31    #kg
global const mₙ     = 1.67492750e-27    #kg
global const MSOLAR = 1.98847e30        #kg

global const LENGTH_TO_KM_CONVERTIONFACTOR = 1e-3
end

module CGSUnits
export G, c, ħ, mₑ, mₙ, MSOLAR, LENGTH_TO_KM_CONVERTIONFACTOR

global const G      = 6.674e-8          #dyne⋅cm²⋅g⁻²
global const c      = 2.99792458e10     #cm⋅s⁻²
global const ħ      = 1.0545919e-27     #erg⋅s
global const mₑ     = 9.10938370e-28    #g
global const mₙ     = 1.67492750e-24    #g
global const MSOLAR = 1.98847e33        #g

global const LENGTH_TO_KM_CONVERTIONFACTOR = 1e-5
end
