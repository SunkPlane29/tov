global const c::Real = G = MSOLAR = 1

global const MASS_UNIT_TO_SI::Real       = 1.989e30
global const LENGTH_UNIT_TO_SI::Real     = 1.477e3
global const TIME_UNIT_TO_SI::Real       = 4.927e-6
global const SI_TO_MASS_UNIT::Real       = MASS_UNIT_TO_SI^(-1)
global const SI_TO_LENGTH_UNIT::Real     = LENGTH_UNIT_TO_SI^(-1)
global const SI_TO_TIME_UNIT::Real       = TIME_UNIT_TO_SI^(-1)

global const PRESSURE_UNIT_TO_SI::Real   = MASS_UNIT_TO_SI*LENGTH_UNIT_TO_SI^(-1)*TIME_UNIT_TO_SI^(-2)
global const SI_TO_PRESSURE_UNIT::Real   = PRESSURE_UNIT_TO_SI^(-1)
global const SI_TO_GEV_FM3::Real         = 6.241509074e-33

global const ħ::Real  = 1.054e-34*SI_TO_MASS_UNIT*SI_TO_LENGTH_UNIT^2*SI_TO_TIME_UNIT^(-1)
global const mₙ::Real = 1.674e-27*SI_TO_MASS_UNIT
global const mₑ::Real = 9.10938370e−31*SI_TO_MASS_UNIT
