
# ---------------------------------------------------------------------------- #
#                     Open Libraries and Define Functions                      #
# ---------------------------------------------------------------------------- #
using SAC, DataFrames, RCall, DSP

# Read the config.txt file to get the
include("bin/read_config.jl")

# ---------------------------------------------------------------------------- #
#                       Set Paths and Define Constants                         #
# ---------------------------------------------------------------------------- #
# Load all of the functions we need
include( function_path*"egf_functions.jl" )
include( function_path*"preprocessing_functions.jl" )
include( function_path*"filtering_functions.jl" )


# Load the list of stations
sta_list = readtable("station_list.txt", header = false, separator = '\t')
lid = sta_list[:,2]

n = size(sta_list, 1)

# ---------------------------------------------------------------------------- #
#                      Get the EGF For Each Station                            #
# ---------------------------------------------------------------------------- #


# Allocate space
xcf = zeros( 2win_len, n)

# Determine the forward and backward operators
bf = lid - ref_sta
ind = (find(lid .== ref_sta) )[1]

# Get the time series for the reference station
ts1 = SAC.read(sta_list[ind,1])

# get the autocorelation for the reference station
if filter_type == "none"
    arr1 = mw_array(ts1.t, win_len, STEP, true)
else
    # define the filter response function
    frf = butt_design(fc, filter_type)
    ts1.t = filt(frf, ts1.t)
end

# Get the cross correlation for the forward, backward and sum of both,
# respective
xcf[:,ind] = arr_xcorr_td(ts1.t, ts1.t)
mlag1[ind] = 2*win_len


find( abs(xcf1[:,ind]) .== maximum(abs(xcf1[:,ind] ) ) )[1]
# Lets cross correlate with the rest of the stations
st_ind = lid[lid .!= ind]
#bf = bf[bf .!= 0]
#st_ind_pos = st_ind[bf .> 0]
#st_ind_neg = st_ind[bf .< 0]

k = 21

if filter_type == "none"
    for i in st_ind
        ts2 = SAC.read(sta_list[i,1])
        xcf[:,i] = arr_xcorr_td( ts1.t, ts2.t )
    end
else
    for i in 1:n
        ts2 = SAC.read(sta_list[i,1])
        xcf[:,i] = arr_xcorr_td( ts1.t, filt(frf, ts2.t) )
    end
end

R"""
    library(fields)
    dev.new()
    image($xcf, col = tim.colors(100) )
"""

#=
xcf_forward = DataFrame(xcf_forward)
xcf_backward = DataFrame(xcf_backward)

writetable("virtual_shot_gather.csv", xcf_forward, separator = ',', header = false)
writetable("virtual_reciever_gather.csv", xcf_backward, separator = ',', header = false)



R"""
    library(fields)
    dev.new()
    image(t($xcf[1:$M,]), col = tim.colors(100) )
"""

=#
