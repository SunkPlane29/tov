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

    return new_eos
end

struct EOS
    pressure::AbstractVector{Real},
    energy_density::AbstractVector{Real},
    eos_interp::Interpolations.Extrapolations,
    eos_fn::Function,
end

function eos_from_file(file::AbstractString, header::AbstractVector ; iscsv::Bool=true)::EOS
    filename = file
    if !iscsv
        filename = TOV.dat2csv(file)
    end

    eos = CSV.File(filename, header=header) |> DataFrame
    eos = remove_negative_slopes(eos)
    #TODO: think about adding other supported units for EoS file
    pressure = eos.p .* MEVFM3_TO_PRESSURE_UNIT
    energy_density = eos.ϵ .* MEVFM3_TO_PRESSURE_UNIT

    eos_interp = linear_interpolation(pressure, energy_density)
    eos_fn(p) = eos_interp(p)

    return EOS(pressure, energy_density, eos_interp, eos_fn)
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
