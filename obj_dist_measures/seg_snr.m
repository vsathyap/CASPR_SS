% function[<output>] = seg_snr(<input>);
%
% Computes per-frame snr as snr_i = 10*log10(||s_i||^2 / ||(s_i-x_i||)
% where s_i is the ith frame of a reference signal s and x_i is the
% ith frame of a reference signals x.
%
% INPUT (in <input> struct):
%
% .s             : Input signal vector (reference signal)
%
% .x             : Processed signal vector (test signal)
%
% .fs (optional) : Sample rate of input signals [Hz]
%                  Default is fs = 20000 Hz.            
%
% .N  (optional) : Length of analysis window [samples]
%                  Default is N = 20 ms
%
% .hop (optional): The number of samples between successive analysis
%                  windows.
%                  Default is hop = 10 ms.
% 
% OUTPUT (in <output> struct):
%
% .ssnr          : Vector with snr per frame.
%
% .tvec          : Vector with sample indeces corresponding to frame
%                  centers.
%
% Author: Jesper Jensen (JSJ), Oticon A/S, 080220.  

function[output] = seg_snr(input);
  
% globals (could be user input as well)
l_limit = -10; %lower limit on perceptually meaningful frame snr
u_limit = 30;  %upper limit on perceptually meaningful frame snr

%default settings
fs_def_Hz = 20000;%sample rate [Hz]
N_def_ms = 20;   %window length in ms.
hop_def_ms = 10; %window hop in ms.

%check validity of inputs
if ~(isfield(input,'s') & isfield(input,'x'))
  disp('error: not enough input arguments.')
  pause;
end

%format signal vectors
s = input.s(:) ; %columnize
x = input.x(:) ; %columnize

if length(s) > length(x)
  disp('warning: truncating signal s')
  s = s(1:length(s));
elseif length(s) < length(x)
  disp('warning: truncating signal x')
  x = x(1:length(x));
end

% fill in for missing inputs
if ~isfield(input,'fs')
  fs = fs_def_Hz;
else
  fs = input.fs;
end

if ~isfield(input,'N')
  N = round(N_def_ms/1000*fs);
else
  N = input.N;
end

if ~isfield(input,'hop')
  hop = round(hop_def_ms/1000*fs);
else
  hop = input.hop;
end

%
% Compute snr per frame
%
NFrames = floor( (length(s) - N)/hop) + 1; % number of frames
ssnr = zeros(NFrames,1);                   % room for results
tvec = zeros(NFrames,1);                   % discete time indeces
evec = zeros(NFrames,1);

rightpointer = N;
for IFrame = 1:NFrames
  index = rightpointer-N+1:rightpointer;   % frame index
  s_frame = s(index);                      % clean frame
  x_frame = x(index);                      % processed frame
  e_frame = x_frame - s_frame;             % error frame

  ssnr(IFrame) = 10*log10(sum(s_frame.^2)/sum(e_frame.^2));
  tvec(IFrame) = round(rightpointer-N/2);  % centre of analysis win
evec(IFrame) = sum(s_frame.^2);

  rightpointer = rightpointer + hop; %move analysis window forward
end

%clip to upper and lower limits for better perceptual relevance
ssnr = min(max(ssnr,l_limit),u_limit);

%export to output struct
output.ssnr = ssnr;
output.tvec = tvec;
evec(find(evec==0)) = eps;
output.evec_dB = 10*log10(evec);
vad_index = find(output.evec_dB>max(output.evec_dB)-30);
output.vad_index = vad_index;