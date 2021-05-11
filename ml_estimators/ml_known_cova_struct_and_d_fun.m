function [lambda_s_ml, lambda_v_ml] = ml_known_cova_struct_and_d_fun(R,Gamma_v,d,ix_ref);
%
% Computes maximum likelihood estimates of lambda_v, lambda_s,
% given d and Gamma_v in the zero-mean Gaussian model:
% Cx = lambda_s*d + lambda_v*Gamma_v
% 
% Equations refer too icassp paper: J. Jensen and M. S. Pedersen,
% "analysis of beamformer directed single-channel noise reduction
% system for hearing aid applications", icassp 2015.
%
% Input:
%        R:       MxM sample covariance matrix of noisy signal
%        Gamma_v: MxM normalized covariance matrix - has "1" in
%                 element in diagonal element (ix_ref,ix_ref).
%        d:       Mx1 relative acoustic transfer function - has "1"
%                 in position  (ix_ref)
%        ix_ref: index of reference microphone
%
% Author: Jesper Jensen, CASPR, AAU, 2019.

d = d/d(ix_ref);%normalize rtf if not already normalized
M = length(R);
iGamma = pinv(Gamma_v);

Qu = eye(M) - 1/(d'*iGamma*d)*d*d'*iGamma;
lambda_v_ml = real(1/(M-1)*trace(Qu*R*iGamma));%eq(3) in UK/JJ
w = iGamma*d/(d'*iGamma*d);%mvdr
lambda_s_ml = real(w'*(R-lambda_v_ml*Gamma_v)*w);%eq(5) in UK/JJ

