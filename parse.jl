function dat2csv(dat_path::AbstractString, csv_path::AbstractString)::AbstractString
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
    ext == ".dat" ||
        throw(ArgumentError("file name doesn't end with `.dat`"))
    return dat2csv(dat_path, "$base.csv")
end

using CSV
using DataFrames

function read_eos_datafile()::DataFrame
    file = "su2njl_eos.csv"
    if !isfile("su2njl_eos.csv")
        dat2csv("su2njl_eos.dat")
    end
    df = CSV.File(file, header=["energy_density", "pressure", "barion_density", "chemical_potential"]) |> DataFrame

    return df
end
