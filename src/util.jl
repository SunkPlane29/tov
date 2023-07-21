using TOV
using DataFrames
using CSV
using Interpolations

function remove_negative_slopes(eos::DataFrame)::DataFrame
    neweos = DataFrame()
    neweos.p = Vector{Float64}()
    neweos.ϵ = Vector{Float64}()

    highest = 0.0
    for (i, p) in enumerate(eos.p)
        if p >= highest
            highest = p
            append!(neweos.p, p)
            append!(neweos.ϵ, eos.ϵ[i])
        end
    end

    return neweos
end

struct EOS
    pressure::AbstractVector{Real}
    energy_density::AbstractVector{Real}
    eos_interp::Interpolations.Extrapolation
    eos_function::Function
end

#NOTE: from this it is clear that the eos file should not already have a header, and, that
#the headers p and ϵ should be in header
function eos_from_file(file::AbstractString, header::AbstractVector{String} ; iscsv::Bool=true)#::EOS
    if !("p" in header) || !("ϵ" in header)
        throw(DomainError("headers p and ϵ must be in header"))
    end

    filename = file
    if !iscsv
        filename = TOV.dat2csv(file)
    end

    eos = CSV.File(filename, header=header) |> DataFrame
    eos = remove_negative_slopes(eos)
    #TODO: think about adding other supported units for EoS file
    pressure = eos.p .* MEVFM3_TO_PRESSURE_UNIT
    energy_density = eos.ϵ .* MEVFM3_TO_PRESSURE_UNIT

    eos_interp = linear_interpolation(pressure, energy_density, extrapolation_bc=Line())
    eos_function(p) = begin
        eos_interp(p)
    end

    return EOS(pressure, energy_density, eos_interp, eos_function)
end

#functions to convert .dat files (commonly from fortran programs) into .csv files, which
#are more machine readable. Stolen from stackoverflow.
function _dat2csv(dat_path::AbstractString, csv_path::AbstractString)::AbstractString
    open(csv_path, write=true) do io
        for line in eachline(dat_path)
            join(io, split(line), ',')
            println(io)
        end
    end
    return csv_path
end

function dat2csv(dat_path::AbstractString)::AbstractString
    base, ext = splitext(dat_path)
    return _dat2csv(dat_path, "$base.csv")
end
