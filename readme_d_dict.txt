Description for loading the dictionary of d-vectors

For exercise 2.i) you are given a dictionary of RATF vectors d_dict.mat. The left hearing aid with the frontal (reference microphone) and rear microphone is used to obtain the RATF vectors.
 
Load the dictionary by including following code in the beginning of the main matlab script:
load(‘d_dict.mat’)
 
An  M x K x J dimensional matrix named ‘d_dict’ will be available in the workspace where M=2 is the number of microphones, K=129 is frequency bin index (STFT window length is 256 and fs = 20 kHz), and J=72 is the angular index. J-index 1 is the RATF vectors associated with the frontal direction (0 degrees), J-index  2 is the RATF vectors associated with direction 5 degree (counter clockwise), J-index 3 is 10 degrees (counter clockwise), and so on.
