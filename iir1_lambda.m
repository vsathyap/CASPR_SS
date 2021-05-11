function lambda = iir1_lambda(tau, fs)
% lambda = iir1_lambda(tau, fs)
%
% Compute coefficient for first-order IIR lowpass filter,
% based on time constant tau [s].
%
% Input:
%  tau: specifies 1st order lowpass filter of the form  
%       b = [lambda,  a= [1 -(1-lambda)]];
%  fs:            Sampling freq.[Hz].
%
% Output:
%  lambda = 1-exp(-1/fs/tau);
%
% 
% Author:
% Jesper Jensen, CASPR, Aalborg University, 2017.
lambda = 1-exp(-(1./fs)./tau);





