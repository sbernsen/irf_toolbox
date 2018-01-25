
# ---------------------------------------------------------------------------- #
#                     Open Libraries and Define Functions                      #
# ---------------------------------------------------------------------------- #
using SAC, DataFrames, RCall, DSP

# Read the config.txt file to get the
config = readtable("config.txt", header = false, separator = ':')

# Assign the values listed in the configuration file
win_len = parse(Int, (config[ config[:,1] .== "window",2])[1] )
STEP = parse(Int, (config[ config[:,1] .== "step", 2] )[1] )
function_path = (config[ config[:,1] .== "function_path", 2])[1]
ref_sta = parse(Int, (config[ config[:,1] .== "reference_station", 2] )[1] )
include(function_path)

# ---------------------------------------------------------------------------- #
#                       Set Paths and Define Constants                         #
# ---------------------------------------------------------------------------- #
# Load the list of stations
sta_list = readtable("station_list.txt", header = false, separator = '\t')
lid = sta_list[:,2]

lsst = lsst[lsst[:,2],:]
n = size(lsst, 1)

# Define the source reciever by location or by filename


# ---------------------------------------------------------------------------- #
#                      Get the EGF For Each Station                            #
# ---------------------------------------------------------------------------- #


# Allocate space
xcf_forward = zeros(win_len, n)
xcf_backward = zeros(win_len, n)

ts_array = zeros(m, n)
forward = true
backward = false



for i = 1:n
    xcfb = zeros(M, 1)
    xcff = zeros(M, 1)
    st_ind = filter!(e->e!=i, int)

    for j in st_ind

        xc = crosscor(ts_array[:,i], ts_array[:,j], [1:(M-1); true] )
        xc = rms_norm(xc)
        xcff = xcff + xc

        xc = crosscor(ts_array[:,j], ts_array[:,i], [1:(M-1); true] )
        xc = rms_norm(xc)
        xcfb = xcfb + xc
    end

    xcf_forward[:,i] = xcff./(n-1)
    xcf_backward[:,i] = xcfb./(n-1)

end

xcf_forward = DataFrame(xcf_forward)
xcf_backward = DataFrame(xcf_backward)

writetable("virtual_shot_gather.csv", xcf_forward, separator = ',', header = false)
writetable("virtual_reciever_gather.csv", xcf_backward, separator = ',', header = false)


#=
R"""
    library(fields)
    dev.new()
    image(t($xcf[1:$M,]), col = tim.colors(100) )
"""
=#
