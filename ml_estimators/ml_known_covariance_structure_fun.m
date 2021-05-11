function [d_ml, lambda_s_ml, lambda_v_ml] = ml_known_covariance_structure_fun(R,Gamma_v,ix_ref);

% Computes maximum likelihood estimates of lambda_v, lambda_s, d,
% given Gamma_v in the zero-mean Gaussian model:
% Cx = lambda_s*d + lambda_v*Gamma_v
%
% Input:
%        R:       MxM sample covariance matrix of noisy signal
%        Gamma_v: MxM normalized covariance matrix - has "1" in
%                 element in diagonal element (ix_ref,ix_ref).
%        ix_ref: index of reference microphone
%
% Author: Jesper Jensen, CASPR, AAU, 2019.


%compute pre-whitened sample covariance matrix
[V,D] = eig(Gamma_v);
Gamma_v_isqrt = V*diag(1./sqrt(diag(D)+eps))*V';
Rtilde = Gamma_v_isqrt*R*Gamma_v_isqrt;
Gamma_v_sqrt = V*sqrt(D)*V';

% EVD 
[VR,DR] = eig(Rtilde);
[eigs_sort,eigs_ix] = sort(real(diag(DR)),'descend');
VR = VR(:,eigs_ix);
DR = diag(eigs_sort);%max eig first

%ML estimation
d_tilde = VR(:,1);
lambda_v_tilde = mean(eigs_sort(2:end));%avg of noise-only eigs
lambda_s_tilde = eigs_sort(1) - lambda_v_tilde;%noisy - noise estimate

% de-whiten -> prepare outputs
d_ml_dewhiten = Gamma_v_sqrt*d_tilde;
d_ml = d_ml_dewhiten/d_ml_dewhiten(ix_ref);%norm to unit ref.
lambda_s_ml = abs(d_ml_dewhiten(ix_ref))^2*lambda_s_tilde;%norm to unit ref
lambda_v_ml = lambda_v_tilde;
