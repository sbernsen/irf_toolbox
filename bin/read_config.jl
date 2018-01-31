#!/usr/bin/env julia

using CSV

# the file is defined on the command line
config = CSV.read(config_file, header = [], delim = ':')

# Assign the values listed in the configuration file
win_len = parse(Int, (config[ config[:,1] .== "window",2])[1] )
STEP = parse(Int, (config[ config[:,1] .== "step", 2] )[1] )
function_path = strip( (config[ config[:,1] .== "function_path", 2])[1] )
ref_sta = parse(Int, (config[ config[:,1] .== "reference_station", 2] )[1] )
filter_type = strip( (config[ config[:,1] .== "filter", 2])[1] )
fc = strip( ( config[ config[:,1].== "corner", 2] )[1] )
fs = parse(Int, (config[ config[:,1] .== "fs", 2] )[1] )

if filter_type == "bp"

    fc = [ parse(Int, split(fc, '/')[1]), parse(Int, split(fc, '/')[2]) ]./fs
elseif filter_type == "none"
    fc = 1
else
    fc = parse(Int, fc[1])
end
