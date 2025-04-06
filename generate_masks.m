close all
clear
clc
%% Select parameters 
% number_of_gaps_in_1sec = 1;
length_of_gap = 5;
[audio_file, fs] = audioread("dataset/IRMAS-Sample/IRMAS-Sample/Training/vio/001__[vio][nod][cou_fol]2194__1.wav");
winLen = 2048; shiftLen = winLen/4; FFTnum = winLen;

% computing necessary info
audio_len_samples = size(audio_file,1);
audio_len_inSec = round(audio_len_samples/fs); % 3 seconds for IRMAS
number_of_frames = ceil(audio_len_samples/shiftLen); % number of hops
frames_perOneSec = ceil(number_of_frames/audio_len_inSec);

% generate mask of ones

mask = ones(winLen/2+1,number_of_frames);

% generate random index within the frame indexes (not on the edge)
%rng('default');
random_index = randi([length_of_gap+2,frames_perOneSec-length_of_gap-2],1,3);

% put gaps of set length in random indexes
for second=1:audio_len_inSec
    starting_index = frames_perOneSec*(second-1)+random_index(second);
    mask(:,starting_index:starting_index+length_of_gap-1) = 0;
end

%% setup DGT based on PHAIN
[win, ~] = generalizedCosWin(winLen, 'hanning');
tight_win = calcCanonicalTightWindow(win, shiftLen);
tight_win = tight_win/norm(tight_win)*sqrt(shiftLen/winLen);
diff_win = numericalDiffWin(tight_win);
zeroPhaseFlag = true;
rotateFlag = true;

[sigIdx, sumIdx, sumArray, ifftArray, rotIdx] = precomputationForFDGT(audio_len_samples, winLen, shiftLen, FFTnum);

%% generate masked spectrogram test
spectrogram_x = FDGT(audio_file, tight_win, sigIdx, FFTnum, rotIdx, zeroPhaseFlag);

% for i=1:size(spectrogram_x,2)
%     spectrogram_x(:,i) = spectrogram_x(:,i)*mask(i);
% end
spectrogram_x = spectrogram_x.*mask;
figure, imagesc(20*log10(abs(spectrogram_x))), colorbar, axis xy
xlabel('Hop number (-)');
ylabel('Frequency (Hz)');
figure, imagesc(angle(spectrogram_x)), colorbar, axis xy

%% Load audio

