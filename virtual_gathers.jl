#!/usr/bin/env julia



# ---------------------------------------------------------------------------- #
#                     Open Libraries and Define Functions                      #
# ---------------------------------------------------------------------------- #
using SAC, DataFrames, RCall, DSP

# Read the config.txt file to get the values of the window length, step size,
# filter configuration, sampling frequency and the reference station

#=
if isdefined(:ARGS)
    config_file = ARGS
else
    error("Input the configuration filepath")
end
=#

config_file = "config.txt"
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
xcf = zeros( 2win_len-1, n)

# Determine the forward and backward operators
bf = lid - ref_sta
ind = (find(lid .== ref_sta) )[1]

# Get the time series for the reference station
ts1 = SAC.read(sta_list[ind,1])

# if filtering is specified, lets define the transfer function
if filter_type != "none"
    # define the filter transfer function
    frf = butt_design(fc, filter_type)
    ts1.t = filt(frf, ts1.t)
end

# get the autocorelation for the reference station
ts1 = mw_array(ts1.t, win_len, STEP, false)
xcf[:, ind] = arr_xcorr_td(ts1, ts1)


# Lets cross correlate with the rest of the stations
st_ind = lid[lid .!= ind]
#bf = bf[bf .!= 0]
#st_ind_pos = st_ind[bf .> 0]
#st_ind_neg = st_ind[bf .< 0]


if filter_type == "none"
    for i in st_ind
        ts2 = SAC.read(sta_list[i,1])
        ts2 = mw_array(ts2, win_len, STEP, false)
        xcf[:,i] = arr_xcorr_td( ts1.t, ts2.t )
    end

else
    for i in st_ind
        ts2 = SAC.read(sta_list[i,1])
        ts2.t = filt(frf, ts2.t)
        ts2 = mw_array(ts2.t, win_len, STEP, false)
        xcf[:,i] = arr_xcorr_td( ts1, ts2 )

#=
        xcf_tmp = arr_xcorr_td( ts1.t, filt(frf, ts2.t) )

        xcf[win_len:(2win_len-1),i] = xcf_tmp[win_len:(2win_len -1)]

        xcf_tmp = arr_xcorr_td( filt(frf, ts2.t), ts1.t )
        xcf[1:win_len,i] = xcf_tmp[1:win_len]
=#
    end

end


xcf = DataFrame(xcf)

writetable("virtual_gather.csv", xcf, separator = ',', header = false)

R"""
source("../plot_virtual_gathers.R")
"""
