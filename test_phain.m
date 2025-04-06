clear
clc
close all

dlen = 44100*10;
addpath(genpath('phase_correction'));
%% parameters
winLen = 2048;
param.w = winLen;
param.a = winLen/4;
param.M = winLen;
param.wtype = 'hann';

%% preparation
data = randn(dlen,1);
param.Ls = length(data);

[win, ~] = generalizedCosWin(param.w, 'hanning');
tight_win = calcCanonicalTightWindow(win, param.a);
tight_win = tight_win/norm(tight_win)*sqrt(param.a/param.w);
zeroPhaseFlag = true;
rotateFlag = true;
[sigIdx, sumIdx, sumArray, ifftArray, rotIdx] = precomputationForFDGT(param.Ls, param.w, param.a, param.M);


%% computation
tic
for cnt=1:10
    spectrogram_z = FDGT(data, tight_win, sigIdx, param.M, rotIdx, zeroPhaseFlag);
end
toc
pause(1)
tic
for cnt=1:10
    signal_z = invFDGT(spectrogram_z, tight_win, sumIdx, sumArray, ifftArray, rotIdx, zeroPhaseFlag)*param.w;
end
toc
