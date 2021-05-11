function [S] = stft(par_stft,in_sigs);

% Input:
%        in_sigs: time domain signals as columns.
%        par_stft: parameters that define the stft
%
% Output:
%        S: numBands x numFrames x M matrix with STFT coefs.
%
%
% Author:
% Jesper Jensen, CASPR, Aalborg University. 2017.

numSamples = size(in_sigs,1);
numFrames = 1 + floor( (numSamples - par_stft.frame_length)/par_stft.D_A);
M = size(in_sigs,2);%number of mics.
numBands = par_stft.N/2+1; %positive freqs.

% make room
S = zeros(numBands,numFrames,M);

for i_m = 1:M
  right_pointer = par_stft.frame_length;
  w_frames = zeros(par_stft.N,numFrames);%room for frames
  for i_frame = 1:numFrames
     index = right_pointer - par_stft.frame_length + 1:right_pointer;%sample index
     w_frames(:,i_frame) = par_stft.awin(:).*in_sigs(index,i_m); %weighted frames
     right_pointer = right_pointer + par_stft.D_A;
  end
  SS =  fft(w_frames);% compute DFT coefs.
  S(:,:,i_m) = SS(1:numBands,:);%save only pos freqs.
end
