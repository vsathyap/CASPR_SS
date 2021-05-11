function win = mod_hann(len);

% Generate Hann window that adds up to one with 50% overlap
% (the standard Hann window does not)
win =(.5 + .5*cos(2*pi*(-(len-1)/2:(len-1)/2)/len))';
