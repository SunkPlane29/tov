#TODO: make util function for opening and interpolating csv file (linear interpolation for
#starters)

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
