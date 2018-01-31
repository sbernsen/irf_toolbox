function [ g, r, m ] = water_deconv( ts1, ts2, wlev )
%water_deconv Computes the impulse response function between two timeseries
%using water level regularization. Refer to Aster et al. (2011)
%'Parameter Estimation and Inverse Problems' for the algorithm.
%
% Input Variables:
%   ts1, ts2 - the m-by-1 timeseries where ts1 is the source and ts2 is the
%              LHS data vector
%   wlev - the water level parameter to compute the deconvolution
% 
% Output Variables:
%   g - the impulse response function in the time domain
%   m - the model norm
%   r - the residual for the respective alpha parameters
%
%Created by Steven Bernsen
%University of Maine
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Pad each timeseries
ts1 = [ ts1; zeros( length(ts1) , 1) ];
ts2 = [ ts2; zeros( length(ts2) , 1) ];

% Compute the fourier transforms 
S = fft(ts1);
X = fft(ts2);


% apply water level damping
for i=1:length(S)
    if abs(S(i)) < wlev
        Sw(i,1)=wlev*(S(i)./abs(S(i)));
    else
        Sw(i,1)=S(i);
    end
end

% get the temporal model and store the necessary information
g=real(ifft(X./Sw));
m=norm(g);
r=norm(real(ifft((X./Sw).*S-X)));

end

