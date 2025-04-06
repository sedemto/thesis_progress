clear
clc
close all

dlen = 44100*10;

%% parameters
winLen = 2048;
param.w = winLen;
param.a = winLen/4;
param.M = winLen;
param.wtype = 'hann';

%% preparation
param.F = frametight(frame('dgtreal', {param.wtype, param.w}, param.a, param.M)); 
%param.F = frametight(frame('dgt', {param.wtype, param.w}, param.a, param.M)); 

data = randn(framelength(param.F, dlen),1);
param.Ls = length(data);
param.F = frameaccel(param.F, param.Ls);


%% computation
tic
for cnt=1:10
    spectrogram_z = framecoef2tf(param.F,frana(param.F, data));
end
toc
pause(1)
tic
for cnt=1:10
    signal_z = frsyn(param.F, frametf2coef(param.F,spectrogram_z));
end
toc
