function [ g, m, r ] = tikh2_deconv( ts1, ts2, alpha)
%tikh0_deconv Computes the impulse response function between two timeseries
%using a second order Tikhonov regularization. Refer to Aster et al. (2011)
%'Parameter Estimation and Inverse Problems' for the algorithm.
%
% Input Variables:
%   ts1, ts2 - the m-by-1 timeseries where ts1 is the source and ts2 is the
%              LHS data vector
%   alpha - the alpha parameter to compute the deconvolution
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
ts1 = [ zeros( floor(length(ts1)/2), 1); ts1; zeros( ceil(length(ts1)/2), 1) ];
ts2 = [ zeros( floor(length(ts2)/2), 1); ts2; zeros( ceil(length(ts2)/2), 1) ]; 

% Compute the fourier transforms 
S = fft(ts1);
X = fft(ts2);



%apply second order Tikhonov regularization
f=fftshift((0:length(S)-1)'/length(S))-0.5;
G=(conj(S).*X)./(conj(S).*S+alpha^2*(f.^4));

% get the temporal model and store the necessary information
g=real(ifft(G));

% To find the best alpha value we will compare the model norm and the
% residual norm
m = norm(G);
r = norm(S.*G-X);

end

