##### -------------------------------------------------------------------- #####
#                              FILTERING_FUNCTIONS                             #
#
# Created by Steven Bernsen because the wheel still isn't good enough
# University of Maine
##### -------------------------------------------------------------------- #####


using DSP

function butt_design( ts, filter_type )
# lowpass_butter applies a low pass, 2 pole butterworth filter to a time series.
# The filter is designed to a normalized sampling frequency of 1 and can be
# scaled afterword. For example, the cutoff frequency of 20 MHz for a
# timeseries sampled at 250 MHz would be fc = 20/250 = 0.08
#
# Input Variables:
#   ts - the input timeseries
#   fc - the corner frequency or frequencies. For bandpass filtering, input as
#       [hp_corner, lp_corner]
#   type - specified as "hp", "lp", "bp" for highpass, lowpass and bandpass,
#       respectively
#
# Output Variables:
#
#
##### -------------------------------------------------------------------- #####

if filter_type == "bp"
    responsetype = Bandpass( fc[1], fc[2] )
elseif filter_type == "lp"
    responsetype = Lowpass(fc)
else
    responsetype = Highpass(fc)
end

designmethod = Butterworth(2)
filter_object = digitalfilter(responsetype, designmethod)


end




#function mtm
#mt_pgram(s; onesided=eltype(s)<:Real, nfft=nextfastfft(n),
#                                fs=1, nw=4, ntapers=iceil(2nw)-1,
#                                window=dpss(length(s), nw, ntapers))
