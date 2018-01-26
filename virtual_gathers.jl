
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

# Determine the forward and backward operators
bf = lid - ref_sta
ind = (find(lid .== 17) )[1]

# Get the time series for the reference station
ts1 = SAC.read(sta_list[ind,1])

# get the autocorelation for the reference station
arr1 = mw_array(ts1.t, win_len, STEP, false)

# Get the cross correlation for the forward, backward and sum of both,
# respective
arr_xc = array_xcorr(arr1, arr1)

# The following two lines are the same for the autocorrelation
xcf_forward[:,ind] = mapslices(mean, arr_xc[1], 2).*ts1.delta # stack the rows
xcf_backward[:,ind] = mapslices(mean, arr_xc[2], 2).*ts1.delta

# Lets cross correlate with the rest of the stations
st_ind = filter!(e->e!=i, ind)

for i in st_ind
    ts2 = SAC.read(sta_list[ind,1])
    arr2 = mw_array(ts2.t, win_len, STEP, false)
    arr_xc = array_xcorr(arr1, arr2)

    xcf_forward[:,i] = mapslices(mean, arr_xc[1], 2).*ts1.delta # stack the rows
    xcf_backward[:,i] = mapslices(mean, arr_xc[2], 2).*ts1.delta

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
