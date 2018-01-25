# ---------------------------------------------------------------------------- #
#                     Open Libraries and Define Functions                      #
# ---------------------------------------------------------------------------- #
using SAC, DataFrames, RCall, DSP


# ---------------------------------------------------------------------------- #
#                       Set Paths and Define Constants                         #
# ---------------------------------------------------------------------------- #
# Load the list of stations
lsst = readtable("station_list.txt", header = false, separator = '\t')
lid = lsst[:,2]

lsst = lsst[lsst[:,2],:]
n = size(lsst, 1)

function_path = "functions/egf_functions.jl"
include(function_path)


# ---------------------------------------------------------------------------- #
#                      Get the EGF For Each Station                            #
# ---------------------------------------------------------------------------- #

# Allocate space
ts = SAC.read(lsst[1,1])
m = length(ts1.t)
M = m-1000
ts_array = zeros(m, n)
xcf_forward = zeros(M, n)
xcf_backward = zeros(M,n)
forward = true
backward = false

int = collect(1:n)

ts_array[:,1] = ts.t
for i = 2:n
        ts = SAC.read(lsst[i,1])
        ts_array[:,i] = ts.t
end

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

R"""
    library(fields)
    dev.new()
    image(t($xcf[1:$M,]), col = tim.colors(100) )
"""
