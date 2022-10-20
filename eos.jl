include("parse.jl")

using DataFrames

struct EOSData
    df::DataFrame
end

# loadEOS() reads the EOS data file containing energy_density, pressure, barion_density
# chemical potential and returns a DataFrame with the data organized in rows and columns
function loadEOS()::EOSData
    df = read_eos_datafile()
    return EOSData(df)
end

# mean() calculates the mean of a given vector
function mean(v::Vector{Real})::Real
    sum = 0.0

    for x in v
        sum += x
    end

    return sum/length(v)
end

# average_slope calculates the mean average of a linear slope from a dataset
function average_slope(eos::EOSData)::Real
    slopes = Vector{Real}()

    pressure_vals = eos.df.pressure
    energydensity_vals = eos.df.energy_density

    for i = 1:length(pressure_vals)-1
        slope = (energydensity_vals[i+1] - energydensity_vals[i])/(pressure_vals[i+1] - pressure_vals[i])
        append!(slopes, slope)
    end

    return mean(slopes)
end

# linearget_energydensity gets the corresponding energy density given the pressure
# assuming the relation is linear. This is a source of error because my curve fitting
# is not ideal and there is some great difference between calculated points and points in the
# dataset
# TODO: look for a better method for this curve fitting (or maybe computing my own EOS)
function linearget_energydensity(avg_slope::Real, ϵ₀::Real, pressure::Real)::Real
    return avg_slope*pressure + ϵ₀
end
