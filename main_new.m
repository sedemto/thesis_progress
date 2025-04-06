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
%soundDir = "dataset/DPAI_originals/audio";
%soundDir = "dataset/sine_wave";
soundDir = "dataset/a08_violin";
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
     
%gaps = 5:5:50 % [ms]
gaps = 1:2; % [1 2 3 4 5] 
M = length(gaps);

N = 2; % # of gaps

methodLabels = {'U_PHAIN'};

for i = 1:length(methodLabels)
    solution.(methodLabels{i}) = cell(NN);  % initialization of restored results
end
SNR  = NaN(NN, M, N, length(methodLabels));  % SNR per gap
TIME = NaN(NN, M, N, length(methodLabels));  % execution time per gap
SNR_procedure = cell(NN, M, N, length(methodLabels));  % iterative behavior

%% parameters

winLen = 2048;
shiftLen = winLen/4;
FFTnum = winLen;

% parameter setting for PHAINs
param.a = shiftLen;
param.M = FFTnum;
param.w = winLen;
param.wtype = 'hann';
param.Ls = length(data{1});
param.g = gabtight(param.wtype,param.a,param.M,param.w);

paramsolver.epsilon = 0.01;  % for stopping criterion

paramsolver.tau = 0.25;  % step size
paramsolver.sigma = 1;  % step size
paramsolver.alpha = 1;  % relaxation parameter

paramsolver.lambda = 1;  % threshold (regularization parameter)

%% test for speed
% signal_x = data{1};
% number = (ceil(size(signal_x,1)/winLen))*winLen;
% times = zeros(1, 100);
% for i=1:100
%     
%     
% 
%     param.F = frametight(frame('dgtreal', {param.wtype, param.w}, param.a, param.M)); 
%     param.F = frameaccel(param.F, param.Ls);
% %     [win, ~] = generalizedCosWin(param.w, 'hanning');
% %     tight_win = calcCanonicalTightWindow(win, param.a);
% %     tight_win = tight_win/norm(tight_win)*sqrt(param.a/param.w);
% %     zeroPhaseFlag = true;
% %     rotateFlag = true;
% %     [sigIdx, sumIdx, sumArray, ifftArray, rotIdx] = precomputationForFDGT(param.Ls, param.w, param.a, param.M);
%     tic;
% %     spectrogram_x= FDGT(signal_x, tight_win, sigIdx, FFTnum, rotIdx, zeroPhaseFlag);
% %     signal_afterIDGT = invFDGT(spectrogram_x, tight_win, sumIdx, sumArray, ifftArray, rotIdx, zeroPhaseFlag)*param.w;
% 
% %     spectrogram_z= dgtreal(signal_x, {'tight', 'hann'}, param.a, param.M);
% %     signal_z = idgtreal(spectrogram_z,{'tight', 'hann'}, param.a, param.M);
%     spectrogram_z = frana(param.F, signal_x);% signal length 128000 samples
%     signal_z = frsyn(param.F, spectrogram_z);
% %     spectrogram_z = reshape(frana(param.F, signal_x),[winLen/2+1,number/shiftLen]);
%     
% 
%     elapsed_time = toc;
%     times(i) = elapsed_time;
% end
% average_time = mean(times)
%% %%%%%% fast DGT %%%%%%%%
[win, ~] = generalizedCosWin(winLen, 'hanning');
tight_win = calcCanonicalTightWindow(win, shiftLen);
tight_win = tight_win/norm(tight_win)*sqrt(shiftLen/winLen);


zeroPhaseFlag = true;
rotateFlag = true;
[sigIdx, sumIdx, sumArray, ifftArray, rotIdx] = precomputationForFDGT(param.Ls, param.w, param.a, param.M);


signal_x = data{1};
% 
% param.G =@(x) FDGT(x, tight_win, sigIdx, FFTnum, rotIdx, zeroPhaseFlag);
spectrogram_x = FDGT(signal_x, tight_win, sigIdx, FFTnum, rotIdx, zeroPhaseFlag);
% signal_tmp1 = invFDGT(spectrogram_x, tight_win, sumIdx, sumArray, ifftArray, rotIdx, zeroPhaseFlag)*winLen;
% spec_fin  =FDGT(signal_tmp1, tight_win, sigIdx, FFTnum, rotIdx, zeroPhaseFlag);
% 
% 
% spectrogram_x(:,40) = 0;
% signal_tmp2 = invFDGT(spectrogram_x, tight_win, sumIdx, sumArray, ifftArray, rotIdx, zeroPhaseFlag)*winLen;
% spec_fin2 = FDGT(signal_tmp2, tight_win, sigIdx, FFTnum, rotIdx, zeroPhaseFlag);
% snr(spec_fin,spec_fin-spec_fin2)
% sum_spec = (spec_fin-spec_fin2)*1000;
% sum_spec(sum_spec.^2 == 0.0) = 0;
% sum_sig = (signal_tmp1-signal_tmp2)*1000;
% a = find(sum_sig);
% %number = (ceil(size(signal_x,1)/winLen))*winLen;
% disp("done")
% spectrogram_z= dgtreal(signal_x, {'tight', 'hann'}, param.a, param.M);
%spectrogram_z= param.F(signal_x);
%spectrogram_z = reshape(frana(param.F, signal_x),[winLen/2+1,number/shiftLen]);
%spectrogram_z = frana(param.F, signal_x);

%plotframe(param.F, spectrogram_z);
% snr(abs(spectrogram_x(:,3:end-2))-abs(spectrogram_z(:,3:size(spectrogram_z,2)-4)),abs(spectrogram_z(:,3:size(spectrogram_z,2)-4)))
%% Printing mag and phase spectrogram of our signal
% figure, imagesc(20*log10(abs(spectrogram_x))), colorbar, axis xy
% xlabel('Hop number (-)');
% ylabel('Frequency (Hz)');
% figure, imagesc(angle(spectrogram_x)), colorbar, axis xy
% 
% figure, imagesc(20*log10(abs(spectrogram_z(:,1:size(spectrogram_z,2)-2)))), colorbar, axis xy
% xlabel('Hop number (-)');
% ylabel('Frequency (Hz)');
% figure, imagesc(angle(spectrogram_z(:,1:size(spectrogram_z,2)-2))), colorbar, axis xy
% signal = invFDGT(spectrogram_x, tight_win, sumIdx, sumArray, ifftArray, rotIdx, zeroPhaseFlag)*winLen;
% 
% spectrogram_y= FDGT(signal, tight_win, sigIdx, FFTnum, rotIdx, zeroPhaseFlag);
% figure, imagesc(20*log10(abs(spectrogram_y))), colorbar, axis xy
% xlabel('Hop number (-)');
% ylabel('Frequency (Hz)');
% figure, imagesc(angle(spectrogram_y)), colorbar, axis xy
% 
% snr(spectrogram_y-spectrogram_x,spectrogram_x)
%% get masks from dpai
% masks = dir("spectrogram_masks\*.mat");
% num_of_masks = length(masks);
% my_masks = zeros(157,num_of_masks);
% for i=1:num_of_masks
%     my_var = "C"+i;
%     tmp = load(masks(i).name).(my_var);
%     my_masks(:,i) = tmp;
% end

%% creating a mask



% create a mask with the size of dgt coef
num_of_coef = size(real(spectrogram_x(1,:)));
my_mask = ones(size(real(spectrogram_x)));

% create a hole in mask
index = num_of_coef(2)-2*round(num_of_coef(2)/3)-3;
my_mask(:,index:index) = 0;

%inpainted signal
gapped_spectrogram = spectrogram_x.*my_mask;
% figure, imagesc(20*log10(abs(gapped_spectrogram))), colorbar, axis xy
% xlabel('Hop number (-)');
% ylabel('Frequency (Hz)');
% figure, imagesc(angle(gapped_spectrogram)), colorbar, axis xy

% D = @(z) z(:,1:end-1) - z(:,2:end);
%spectrogram_diff = FDGT(gapped_spectrogram, diff_win, sigIdx, param.M, rotIdx, zeroPhaseFlag);
%omega = calcInstFreq(gapped_spectrogram, spectrogram_diff, param.M, param.w, rotateFlag);
%gapped_spectrogram = instPhaseCorrection(gapped_spectrogram, omega, param.a, param.M);

% figure, imagesc(20*log10(abs(gapped_spectrogram))), colorbar, axis xy
% figure, imagesc(angle(gapped_spectrogram)), colorbar, axis xy
% 
% gapped_signal = invFDGT(gapped_spectrogram, tight_win, sumIdx, sumArray, ifftArray, rotIdx, zeroPhaseFlag)*winLen;
% figure, plot(gapped_signal)
% xlabel('Samples (-)');
% ylabel('Amplitude');
% gapped = FDGT(gapped_signal, tight_win, sigIdx, param.M, rotIdx, zeroPhaseFlag);
% figure, imagesc(20*log10(abs(gapped))), colorbar, axis xy
% xlabel('Hop number (-)');
% ylabel('Frequency (Hz)');
% figure, imagesc(angle(gapped)), colorbar, axis xy
% %audiowrite("afterDGTinvDGT.wav",gapped_signal, fs);
% snr(gapped_spectrogram(:,1:100)-gapped(:,1:100),gapped_spectrogram(:,1:100))
%% testing
% assuming all signals are equal length
% for nn=1:NN
% 
%     current_signal = data{nn};
%     spect = param.G(current_signal);
%     
%     for m=1:size(my_masks,2)
% 
%         % replace masks if incompatible legth with spec
%         if size(spect,2)~= size(my_masks(:,m),1)
%             my_masks = my_masks(1:size(spect,2),:);
%             disp("reduced size")
%         end
% 
%         current_mask  = my_masks(:,m)';
%         gapped_spec = spect.*current_mask;
% 
% %         figure, imagesc(20*log10(abs(gapped_spec))), colorbar, axis xy
% %         figure, imagesc(angle(gapped_spec)), colorbar, axis xy
%         
%         % save gapped spectrogram as solution
%         solution.(methodLabels{1}){nn,m} = gapped_spec;
%         
%         % next part depends on overlap of the windows (for our case 75 %)
%         % that means we need to calculate all windows affected by gap
%         % do it for all gaps
%         indxs_gaps = find(all(gapped_spec == 0,1));
%         
%         for gap=1:5
%             gap = gap*m;
%             % get index of gap and also affected windows
%             gap_indx = indxs_gaps(gap-m+1:gap);
%             segment_indx = (gap_indx(1)-3:gap_indx(1)+m-1+3);
% 
%             segment.mask = current_mask(segment_indx);
%             segment.data = spect(:,segment_indx);
%             segment.gapped = gapped_spec(:,segment_indx);
%             segment.max = max(abs(segment.data));
%         end
%         %segments = 
%         for i=1:length(indxs_gaps)
% 
%         end
%         
%         
%     end
%     
%     
%     
% 
% 
% end
%% inpainting

indeces = [index];
full.length =  length(spectrogram_x(1,:));
segment.length = 8;
for i = 1:length(methodLabels)
    solution.(methodLabels{i}){1} = gapped_spectrogram;
end
for n = 1:N
    indxs_gap = [indeces(n):indeces(n)+1];
    if length(indxs_gap) ~= segment.length
        to_pad = segment.length - length(indxs_gap);
        to_pad_each = 4;
        indxs_gap = [(indeces(n)-to_pad_each):(indeces(n)+1+to_pad_each)];
    end
    disp(indxs_gap)
    fprintf('\nGap Number: %d / %d\n', n, N)
    % making a gap
    segment.mask = my_mask(:,indxs_gap);
    segment.data = spectrogram_x(:,indxs_gap);
    segment.max = max(abs(segment.data));
%     segment.gapped = segment.data.*segment.mask;
    [sigIdx2, sumIdx2, sumArray2, ifftArray2, rotIdx2] = precomputationForFDGT(size(segment.data,2)*shiftLen, winLen, shiftLen, FFTnum);
%     temp = invFDGT(segment.data, tight_win, sumIdx2, sumArray2, ifftArray2, rotIdx2, zeroPhaseFlag)*winLen;
%     figure, plot(temp)
%     temp_max = max(abs(temp));
%     disp(temp_max)
%     new_data = temp/temp_max;
%     figure, plot(new_data)
%     
%     segment.data = FDGT(new_data, tight_win, sigIdx2, FFTnum, rotIdx2, zeroPhaseFlag);

    segment.gapped = segment.data.*segment.mask;

%     figure, imagesc(20*log10(abs(segment.gapped))), colorbar, axis xy
%     figure, imagesc(angle(segment.gapped)), colorbar, axis xy
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% U-PHAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf('U-PHAIN...\n')
    tic
    param.type = 'U';
    paramsolver.I = 100;
    paramsolver.J = 10;
    
    [segment.solution, SNR_procedure{nn, 1, n, 1}] = UPHAIN_ltfat(segment.gapped, segment.mask, param, paramsolver, segment.data);
    
    data = invFDGT(segment.data, tight_win, sumIdx2, sumArray2, ifftArray2, rotIdx2, zeroPhaseFlag)*winLen;
    data_gapped =invFDGT(segment.gapped, tight_win, sumIdx2, sumArray2, ifftArray2, rotIdx2, zeroPhaseFlag)*winLen;
%     figure;plot(data,'Color', [0 0 1 1]); hold on
%     plot(data_gapped,'Color', [0 1 1 0.5])
    solution.U_PHAIN{1}(:,indxs_gap) = segment.solution;
    result = invFDGT(solution.U_PHAIN{1}(:,indxs_gap), tight_win, sumIdx2, sumArray2, ifftArray2, rotIdx2, zeroPhaseFlag)*winLen;

%     plot(result,'Color', [1 0 0 0.2]);
    snr(segment.data, segment.data-solution.U_PHAIN{1}(:,indxs_gap))
    snr(data, data-data_gapped)
    snr(data,data-result)
    hold off
    %alpha(a,0.1)
    TIME(nn, m, n, 6) = toc;

end
%first = invFDGT(segment.gapped, tight_win, sumIdx, sumArray, ifftArray, rotIdx, zeroPhaseFlag)*winLen;
%segment.solution = FDGT(first, tight_win, sigIdx, FFTnum, rotIdx, zeroPhaseFlag);
%figure, imagesc(20*log10(abs(segment.solution))), colorbar, axis xy

%% 
figure, imagesc(20*log10(abs(solution.U_PHAIN{1}))), colorbar, axis xy
figure, imagesc(angle(solution.U_PHAIN{1})), colorbar, axis xy
ending = invFDGT(solution.U_PHAIN{1}, tight_win, sumIdx, sumArray, ifftArray, rotIdx, zeroPhaseFlag)*winLen;
%audiowrite("output_sine_phainTF.wav", ending,fs)
gap1 = round(length(signal_x)/size(solution.U_PHAIN{1},2).*(indeces(2)+15));
gap2 = round(length(signal_x)/size(solution.U_PHAIN{1},2).*(indeces(2)+20));
snr(gapped_signal(5000:6000), gapped_signal(5000:6000)-ending(5000:6000))
calcSNR(signal_x(gap1:gap2), ending(gap1:gap2))
%% plot

snr_vec = squeeze(median(SNR, [1,3]));
%figure(Position = [614 157 873 830])
p = plot(snr_vec, LineWidth = 2);
grid on
legend(["U-PHAIN"], Interpreter = 'latex')
xlabel('gap length 5:5:50 [ms]', Interpreter = 'latex')
xticks = 1:10;
xticklabels = ["5", "10", "15", "20", "25", "30", "35", "40", "45", "50"];
ylabel('SNR at gaps [dB]', Interpreter = 'latex')
ax = gca;
ax.FontSize = 15;
ax.TickLabelInterpreter = 'latex';
