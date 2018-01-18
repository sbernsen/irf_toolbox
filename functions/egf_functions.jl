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
