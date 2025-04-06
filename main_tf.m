%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%         PHASE-AWARE AUDIO INPAINTING in TF domain (PHAIN-TF)            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc

rng(0)

%% paths

addpath(genpath('dataset'))
addpath('utils')
addpath(genpath('phase_correction'));
addpath('PHAIN');
%% loading

%
%soundDir = "dataset/IRMAS_five_seconds/audio_original_example00";
soundDir = "dataset/DPAI_originals/audio";
%soundDir = "dataset/sine_wave";
%soundDir = "dataset/a08_violin";
ext = ".wav";
Sounds = dir(soundDir + "*" + ext);
NN = length(Sounds);
data = cell(NN,1);
info = audioinfo(Sounds(1).name);
fs = info.SampleRate;
for nn = 1:NN
    data{nn} = audioread(Sounds(nn).name);
end
clear audio info



%% settings
methodLabels = {'U_PHAIN'};
for i = 1:length(methodLabels)
    solution.(methodLabels{i}) = cell(NN);  % initialization of restored results
end
SNR_spec  = NaN(NN, 6, length(methodLabels)); % 6 masks, NN=number of sigs
SNR_sig  = NaN(NN, 6, length(methodLabels));
TIME = NaN(NN, 6, length(methodLabels));
RESULTS =NaN(NN, 6, length(methodLabels));
%% parameters

winLen = 2048; % window length
shiftLen = winLen/4; % hop size
FFTnum = winLen; % number of rows

% parameter setting for PHAINs
param.a = shiftLen;
param.M = FFTnum;
param.w = winLen;
param.wtype = 'hann'; % window type
param.Ls = length(data{1}); % length of input audio signal

paramsolver.epsilon = 0.01;  % for stopping criterion
paramsolver.tau = 0.25;  % step size
paramsolver.sigma = 1;  % step size
paramsolver.alpha = 1;  % relaxation parameter
paramsolver.lambda = 0.01;  % threshold (regularization parameter)

%% fast DGT and invDGT
% creating a tight window
[win, ~] = generalizedCosWin(winLen, param.wtype);
tight_win = calcCanonicalTightWindow(win, shiftLen);
tight_win = tight_win/norm(tight_win)*sqrt(shiftLen/winLen);

% precomputation for fast DGT
zeroPhaseFlag = true;
[sigIdx, sumIdx, sumArray, ifftArray, rotIdx] = precomputationForFDGT(param.Ls, param.w, param.a, param.M);

% fast DGT and invDGT definition
param.G =@(x) FDGT(x, tight_win, sigIdx, FFTnum, rotIdx, zeroPhaseFlag);
param.G_adj =@(u) invFDGT(u, tight_win, sumIdx, sumArray, ifftArray, rotIdx, zeroPhaseFlag)*winLen;
%% get masks from dpai
masks = dir("spectrogram_masks\*.mat");
num_of_masks = length(masks);
my_masks = zeros(157,num_of_masks);
for i=1:num_of_masks
    my_var = "C"+i;
    tmp = load("spectrogram_masks\"+masks(i).name).(my_var);
    my_masks(:,i) = tmp;
end
%% testing
% assuming all signals are equal length
for nn=1:NN % iterate signals

    current_signal = data{nn};
    spect = param.G(current_signal);
    
    for m=1:size(my_masks,2) % iterate masks 1--6 from janssen2

        % replace masks if incompatible length with spec
        if size(spect,2)~= size(my_masks(:,m),1)
            my_masks = my_masks(1:size(spect,2),:);
            disp("reduced size")
        end

        current_mask  = my_masks(:,m)';
        gapped_spec = spect.*current_mask;

%         figure, imagesc(20*log10(abs(gapped_spec))), colorbar, axis xy
%         figure, imagesc(angle(gapped_spec)), colorbar, axis xy
        
        % save gapped spectrogram as solution
        solution.(methodLabels{1}){nn,m} = gapped_spec;
        
        % next parts depends on overlap of the windows (for our case 75 %)
        % that means we need to calculate all windows affected by gap
        % do it for all gaps 
        indxs_gaps = find(all(gapped_spec == 0,1));
        tic
        for gap=1:5 % iterate gaps
            gap_ = gap*m;
            % get index of gap and also affected windows
            gap_indx = indxs_gaps(gap_-m+1:gap_);
            %disp(gap_indx(end))
            pad = 4; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% needs proper gap restriction, 12 for testing
            pad = fix_pad(pad,gap_indx);
            segment_indx = (gap_indx(1)-pad:gap_indx(end)+pad);

            segment.mask = current_mask(segment_indx);
            segment.data = spect(:,segment_indx);
%%%%%
%             [sigIdx2, sumIdx2, sumArray2, ifftArray2, rotIdx2] = precomputationForFDGT(size(segment.data,2)*shiftLen, winLen, shiftLen, FFTnum);
%             temp = invFDGT(segment.data, tight_win, sumIdx2, sumArray2, ifftArray2, rotIdx2, zeroPhaseFlag)*winLen;
%             %figure, plot(temp)
%             temp_max = max(abs(temp));
%             %disp(temp_max)
%             new_data = temp/temp_max;
%             %figure, plot(new_data)
%              
%             segment.data = FDGT(new_data, tight_win, sigIdx2, FFTnum, rotIdx2, zeroPhaseFlag);
 %%%%
            segment.gapped = segment.data.*segment.mask;%gapped_spec(:,segment_indx);
            fprintf('U-PHAIN...gap number: %d\n',gap)
            
            param.type = 'U';
            paramsolver.I =100;
            paramsolver.J = 10;
            [segment.solution] = PHAINmain_TF(segment.gapped, segment.mask, param, paramsolver,segment.data);
%             figure, imagesc(20*log10(abs(segment.data))), colorbar, axis xy
%             figure, imagesc(20*log10(abs(segment.gapped))), colorbar, axis xy
%             figure, imagesc(20*log10(abs(segment.solution))), colorbar, axis xy
%             temp_sig = invFDGT(segment.solution, tight_win, sumIdx2, sumArray2, ifftArray2, rotIdx2, zeroPhaseFlag)*winLen*temp_max;
%             segment.solution = FDGT(temp_sig, tight_win, sigIdx2, FFTnum, rotIdx2, zeroPhaseFlag);
            solution.(methodLabels{1}){nn,m}(:,segment_indx) = segment.solution;

        end
        TIME(nn,m,1) = toc;
        SNR_spec(nn,m,1) = snr(spect,spect-solution.(methodLabels{1}){nn,m});
%         figure, imagesc(20*log10(abs(spect))), colorbar, axis xy
%         figure, imagesc(angle(spect)), colorbar, axis xy
% 
%         figure, imagesc(20*log10(abs(gapped_spec))), colorbar, axis xy
%         figure, imagesc(angle(gapped_spec)), colorbar, axis xy
% 
%         figure, imagesc(20*log10(abs(solution.(methodLabels{1}){nn,m}))), colorbar, axis xy
%         figure, imagesc(angle(solution.(methodLabels{1}){nn,m})), colorbar, axis xy
        fprintf('done for example: %d, mask: %d!\n',nn,m)
        signal_data = param.G_adj(spect);
        signal_result = param.G_adj(solution.(methodLabels{1}){nn,m});
%         figure;plot(signal_data,'Color',[0 0 1 1]), hold on
%         plot(signal_result,'Color',[1 0 0 0.5])
        name_audio = "results/example"+nn+"mask"+m+".wav";
        %RESULTS(nn,m,1) = 
        audiowrite(name_audio,signal_result,fs);
        SNR_sig(nn,m,1) = snr(signal_data,signal_data-signal_result);
    end
end