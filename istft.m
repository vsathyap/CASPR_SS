function [s] = istft(par_stft,S_stft);

% Input:
%        S_stft: numBands x numFrames matrix with STFT coefs.
%                (numBands = N/2+1, where N is fft order)
%        par_stft: parameters that define the stft
%
% Output:
%        s: time-domain signal.
%
% Author:
%        Jesper Jensen, CASPR, Aalborg University, 2017.


[numBands,numFrames] = size(S_stft);
numSamples = par_stft.N+par_stft.D_A*(numFrames - 1);
s = zeros(numSamples,1);

%expand with negative freqs.
S_stft_tot = zeros(par_stft.N,numFrames);
S_stft_tot(1:numBands,1:numFrames) = S_stft;
S_stft_tot(numBands+1:par_stft.N,1:numFrames) = flipud(conj(S_stft(2:end-1,:)));
s_frames = ifft(S_stft_tot);%compute reconstructed frames

% now overlap-add
right_pointer = par_stft.N;
for iFrame = 1:numFrames
    index = right_pointer - par_stft.frame_length + 1:right_pointer;%sample index
    
    s(index) = s(index) + par_stft.swin(:).*s_frames(:,iFrame); %weighted frames
    right_pointer = right_pointer + par_stft.D_A;
end
