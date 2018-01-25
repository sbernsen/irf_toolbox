using DSP


function padded_mw_array(ts, win_len, STEP)
# Assemble an array for the tukey window where each column is a new step in the
# time series then pad the array with zeros.
#
# Input Variables:
#   ts - the m-by-1 time series
#   win_len - the integer value of the window length
#   STEP - the integer step size for the number of
#
# Output Variables:
#
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Set the tukey alpha value; 0 - rectangle window, 1 - Hann window
alpha = 0.2

# create the tukey window
win = tukey(win_len, alpha)

# Lets get the integer values of the time series for each window
a = [ 1 ]
b = [ win_len ]
m = length(ts)
i = 1


while (b[i]+STEP) <= m
    a = push!(a, a[i] + STEP )
    b = push!(b, b[i] + STEP )

    i = i + 1
end

n = length(a)
padded_array = zeros(2*win_len, n)

for i in 1:n
    padded_array[1:win_len, i] = ts[ a[i]:b[i] ].*win
end

return padded_array

end


function array_xcorr(arr1, arr2)
# Compute the cross correlation in the frequency domain of each column
#
# Input Variables:
#   arr1, arr2 - padded arrays of equal dimensions to be cross correlated
#
# Output Variables:
#   arr_xc - the cross correlation of equal dimensions the input arrays
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Compute the FFT on each column
arr1 = mapslices(fft, arr1, 1)
arr2 = mapslices(fft, arr2, 1)

# Compute the dot product
fft_arr = arr1.*conj(arr2)

# Get the inverse FFT
arr_xc = real( mapslices( ifft, fft_arr, 1) )

return arr_xc, fft_arr

end

function xcorr_td(ts1, ts2, forward)
# compute the time domain cross correlation for two equal length time series for positive or negative lags
#
# Input Variables:
#	ts1, ts2 - the m-by-1 vectors to cross correlate
#	forward - boolean value (true/false) to compute the xcf in positive or negative lags
# Output Variables:
#	xcf - the 2m-by-1 cross correlation vector

# Allocate space
n = length(ts1)
xcf = zeros(n-1, 1)

# Compute the cross correlation
if forward == "true"
	for i = 1:n-1
		xcf[i] = ( (ts1[1:(n-i+1)]'*ts2[i:n])/(std(ts1[1:(n-i+1)])*std(ts2[i:n]) ) )[1]
	end
else
	for i = 1:n-1
		xcf[i] = (  (ts2[1:(n-i+1)]'*ts1[i:n])/(std(ts2[1:(n-i+1)])*std(ts1[i:n]) ) )[1]
	end
end

return xcf

end


function rms_norm(ts)
# Compute the cross correlation in the frequency domain of each column
#
# Input Variables:
#
# Output Variables:
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

n = length(ts)
rms = sqrt( (1/n)*(ts'*ts) )

ts = ts./rms

return ts

end


################################################################################
################################# IO Functions #################################
################################################################################


function matrix2SAC(data_arr, sta_list)
#
#
# This function was made specific for 'virtual_gathers.jl' because each gather
# corresponds to the list of stations used to create them.
end
