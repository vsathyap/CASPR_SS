% function d_vecs = compute_d_vecs(par,look_dir_deg,i_mics,i_ref)
%
% Computes relative transfer functions wrt. reference microphone
% for specified target direction. This is done by estimating
% the cross-power spectral density matrix Cs for each freq.
% band for a ssn point target in the desired direction, 
% and using the eigen vector corresponding to the largest
% eigen value as an estimate of the transfer function.
%
% Input:
%        pa: struct with setting parameters  
%        look_dir_deg: desired look direction in [deg] (0 is frontal, counting
%                      anti-clockwise)
%        i_mics: index of mics used in this array. Microphone nomenclature
%                -  1: left, front, 2: left, rear, 3: right,front, 4: right, rear. 
%        i_ref: index of reference mic.
% 
% Output:
%        d_vecs: M x numBands matrix with relative acoustic transfer
%                function vectors as columns.
%
% Author:
%        Jesper Jensen, CASPR, Aalborg University, 2017.

function d_vecs = compute_d_vecs(par,look_dir_deg,i_mics,i_ref)
%
% to find d-vectors, use ssn for calibration.
%

% pre-stored file with mic.sigs for ssn point noise sources for different
% angles.
mic_sig_filename = 'mic_sigs_ssn.mat';

% load microphone signals for all directions
load_cmd = ['load ' par.sim.micsigDir mic_sig_filename]; 
eval(load_cmd) %load Xo_ssn and corresponding theta

% pick out signals from target direction.
look_dir_vec = theta/2/pi*360;%vector theta loaded just above
%i_theta = find(look_dir_vec == look_dir_deg);%find index of target direction
[val,i_theta] = min(abs(look_dir_vec - look_dir_deg));

% stack all mic sigs for this direction
all_mic_sig_mat = [Xo_ssn.frontLeft.x(:,i_theta) Xo_ssn.rearLeft.x(:,i_theta) Xo_ssn.frontRight.x(:,i_theta) Xo_ssn.rearRight.x(:,i_theta)];

% get mic sigs for specified mic.array
mic_sig_mat = all_mic_sig_mat(:,i_mics);

% compute STFT of all mic.sigs.
S = stft(par.stft,mic_sig_mat);
[numBands, numFrames,M] = size(S);

% compute d_vecs
ii_ref = find(i_mics == i_ref);
d_vecs = zeros(M,numBands);%room for output

for i_band = 1:numBands
% compute covariance matrix
     data_matrix = squeeze(S(i_band,:,:)).';
     Cs = 1/numFrames*data_matrix*data_matrix';

     % find look vector
     [VV,DD] = eig(Cs);
     [y,i]=max(diag(DD));%compute d as eigen vector for max eigen value
     d = VV(:,i);
     d_vecs(:,i_band) = d/d(ii_ref);%relative transfer function wrt. ref. mic
end
