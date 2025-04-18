%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                   PHASE-AWARE AUDIO INPAINTING (PHAIN)                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear
clc

rng(0)

%% paths

addpath(genpath('dataset'))
addpath('utils')
addpath(genpath('phase_correction'));
addpath('PHAIN');

%% loading
%soundDir = "dataset/a08_violin";
%soundDir = "dataset/a16_clarinet";
%soundDir = "dataset/sine_wave";
soundDir = "dataset/DPAI_originals/audio_original_example6";
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
     
gaps = 20:5:30; % [ms]
M = length(gaps);

N = 3; % # of gaps

%methodLabels = {'B_PHAIN', 'B_PHAIN_oracle', 'R_PHAIN', 'R_PHAIN_oracle', 'UR_PHAIN', 'U_PHAIN'};
methodLabels = {'U_PHAIN'};
for i = 1:length(methodLabels)
    solution.(methodLabels{i}) = cell(NN, M);  % initialization of restored results
end

SNR  = NaN(NN, M, N, length(methodLabels));  % SNR per gap
TIME = NaN(NN, M, N, length(methodLabels));  % execution time per gap
SNR_procedure = cell(NN, M, N, length(methodLabels));  % iterative behavior

%% parameters

winLen = 2800;
shiftLen = winLen/4;
FFTnum = winLen;

% parameter setting for PHAINs
param.a = shiftLen;
param.M = FFTnum;
param.w = winLen;
param.wtype = 'hann';
%param.g = gabtight(param.wtype,param.a,param.M,param.w);
%param.F = frametight(frame('dgtreal', {param.wtype, param.w}, param.a, param.M));
paramsolver.epsilon = 0.01;  % for stopping criterion

paramsolver.tau = 0.25;  % step size
paramsolver.sigma = 1;  % step size
paramsolver.alpha = 1;  % relaxation parameter

paramsolver.lambda = 1;  % threshold (regularization parameter)

%% inpainting

for nn = 1:NN

    signal = data{nn};

    for m = 1:M

        fprintf('\nSignal: %d / %d', nn, NN)
        fprintf('\nGap Length: %d [ms]\n', gaps(m))

        gapLength = gaps(m); % [ms]
        h = round(fs*gapLength/1000); % [samples]
        full.length = length(signal);
        full.mask = true(full.length, 1);

        notDegraded = 0.5; % [s]
        segment.length = round((length(signal) - 2*notDegraded*fs) / N);

        for i = 1:length(methodLabels)
            solution.(methodLabels{i}){nn, m} = signal;
        end

        for n = 1:N
            
            fprintf('\nGap Number: %d / %d\n', n, N)
            idx = round(notDegraded*fs) + ((n - 1)*segment.length+1:n*segment.length);
%%PP: cast making a gap prerabat urcite, lebo sa to riesi v case
            % making a gap      
            s = round((winLen + 1) + rand()*(segment.length - 2*winLen - h));
            f = s + h - 1;
            segment.mask = true(segment.length, 1);
            segment.mask(s:f) = false;
            full.mask(idx) = segment.mask;
            
            segment.data_orig = signal(idx);
            
            segment.max = max(abs(segment.data_orig));
            segment.data = segment.data_orig;%/segment.max;
            segment.gapped = segment.data.*segment.mask;
            segment.data_max = segment.data_orig/segment.max;
            segment.gapped_max = segment.data_max.*segment.mask;

            % shortening the segment
            [firstIdx, L] = shortenForDGT(winLen, shiftLen, s, f);
            origL = L;
            enoughL = ceil(L/lcm(shiftLen, FFTnum))*lcm(shiftLen, FFTnum);
            if L < enoughL
                L = enoughL;
            end

            lastIdx = firstIdx + L - 1;
            if firstIdx >= 1 && lastIdx <= segment.length
                segment.mask = segment.mask(firstIdx:lastIdx);
                segment.gapped = segment.gapped(firstIdx:lastIdx);
                segment.data = segment.data(firstIdx:lastIdx);

                segment.gapped_max = segment.gapped_max(firstIdx:lastIdx);
                segment.data_max = segment.data_max(firstIdx:lastIdx);
                idx = idx(firstIdx:firstIdx + L - 1);
                st = 1;
                ed = L;
            else
                firstIdx = max(firstIdx, 1);
                padding = zeros(lastIdx - length(segment.data), 1);
                segment.mask = [segment.mask(firstIdx:end); true(size(padding))];
                segment.gapped = [segment.gapped(firstIdx:end); padding];
                segment.data = [segment.data(firstIdx:end); padding];

                segment.gapped_max = [segment.gapped_max(firstIdx:end); padding];
                segment.data_max = [segment.data_max(firstIdx:end); padding];
                idx = idx(firstIdx:firstIdx + length(segment.data) - length(padding) - 1);
                st = 1;
                ed = L - length(padding);
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%% B-PHAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%

%             fprintf('B-PHAIN...\n')
%             tic
%             param.type = 'B';
%             paramsolver.I = 1000;
% 
%             [segment.solution, SNR_procedure{nn, m, n, 1}] = ...
%                 PHAINmain(segment.gapped, segment.mask, param, paramsolver, segment.data);
% 
%             solution.B_PHAIN{nn, m}(idx) = segment.solution(st:ed)*segment.max;
%             TIME(nn, m, n, 1) = toc;
% 
% 
%             %%%%%%%%%%%%%%%%%%%%%%% B-PHAIN (oracle) %%%%%%%%%%%%%%%%%%%%%%
% 
%             fprintf('B-PHAIN (oracle)...\n')
%             tic
%             param.type = 'Bora';
%             paramsolver.I = 1000;
% 
%             [segment.solution, SNR_procedure{nn, m, n, 2}] = ...
%                 PHAINmain(segment.gapped, segment.mask, param, paramsolver, segment.data);
% 
%             solution.B_PHAIN_oracle{nn, m}(idx) = segment.solution(st:ed)*segment.max;
%             TIME(nn, m, n, 2) = toc;
% 
% 
%             %%%%%%%%%%%%%%%%%%%%%%%%%%% R-PHAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%             fprintf('R-PHAIN...\n')
%             tic
%             param.type = 'R';
%             paramsolver.I = 100;
%             paramsolver.J = 10;
% 
%             [segment.solution, SNR_procedure{nn, m, n, 3}] = ...
%                 PHAINmain(segment.gapped, segment.mask, param, paramsolver, segment.data);
% 
%             solution.R_PHAIN{nn, m}(idx) = segment.solution(st:ed)*segment.max;
%             TIME(nn, m, n, 3) = toc;
% 
% 
%             %%%%%%%%%%%%%%%%%%%%%%% R-PHAIN (oracle) %%%%%%%%%%%%%%%%%%%%%%
% 
%             fprintf('R-PHAIN (oracle)...\n')
%             tic
%             param.type = 'Rora';
%             paramsolver.I = 100;
%             paramsolver.J = 10;
% 
%             [segment.solution, SNR_procedure{nn, m, n, 4}] = ...
%                 PHAINmain(segment.gapped, segment.mask, param, paramsolver, segment.data);
% 
%             solution.R_PHAIN_oracle{nn, m}(idx) = segment.solution(st:ed)*segment.max;
%             TIME(nn, m, n, 4) = toc;
% 
% 
%             %%%%%%%%%%%%%%%%%%%%%%%%%% UR-PHAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%             fprintf('UR-PHAIN...\n')
%             tic
%             param.type = 'UR';
%             paramsolver.I = 100;
%             paramsolver.J = 10;
% 
%             [segment.solution, SNR_procedure{nn, m, n, 5}] = ...
%                 PHAINmain(segment.gapped, segment.mask, param, paramsolver, segment.data);
% 
%             solution.UR_PHAIN{nn, m}(idx) = segment.solution(st:ed)*segment.max;
%             TIME(nn, m, n, 5) = toc;
% 

            %%%%%%%%%%%%%%%%%%%%%%%%%%% U-PHAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fprintf('U-PHAIN...\n')
            tic
            param.type = 'U';
            paramsolver.I = 100;
            paramsolver.J = 10;
            %figure,plot(segment.data);
            %figure,plot(segment.gapped);
            [segment.solution, SNR_procedure{nn, m, n, 1}] = ...
                PHAINmain(segment.gapped, segment.mask, param, paramsolver, segment.data);
            
            [segment.solution_max, SNR_procedure{nn, m, n, 1}] = ...
                PHAINmain(segment.gapped_max, segment.mask, param, paramsolver, segment.data_max);
                    
            solution.U_PHAIN{nn, m}(idx) = segment.solution(st:ed);%*segment.max;
            solution2 = segment.solution_max(st:ed)*segment.max;
            solution1 = segment.solution(st:ed);
            inv_mask =  (1-segment.mask);
            indx_gap  = find(inv_mask);

            figure, plot(segment.data(indx_gap)); hold on
            plot(solution1(indx_gap),'Color',[0 0 1 1]);
            plot(solution2(indx_gap),'Color',[1 0 0 0.5]);
%             figure, plot(solution.U_PHAIN{nn, m}(idx));
            snr(segment.data(indx_gap), segment.data(indx_gap)-solution1(indx_gap))
            snr(segment.data(indx_gap), segment.data(indx_gap)-solution2(indx_gap))
%             snr(abs(signal(idx)), abs(signal(idx)-solution.U_PHAIN{nn, m}(idx)))
            

            TIME(nn, m, n, 1) = toc;

        end
   
        % calculating SNR
        fprintf('\nevaluation start\n')
        for i = 1:length(methodLabels)
            restored = solution.(methodLabels{i}){nn, m};
            groundtruth = signal(~full.mask);
            result = restored(~full.mask);
            for n = 1:N
                SNR(nn, m, n, i) = calcSNR(groundtruth(1 + h*(n - 1):h*n), result(1 + h*(n - 1):h*n));
                snr(signal(idx),signal(idx)-solution.U_PHAIN{nn, m}(idx))
            end
        end
        fprintf('\nevaluation done!\n')

    end

end

%% plot

snr_vec = squeeze(median(SNR, [1,3]));
figure(Position = [614 157 873 830])
p = plot(snr_vec, LineWidth = 2);
grid on
legend(["U-PHAIN"], Interpreter = 'latex')
xlabel('gap length 5:5:50 [ms]', Interpreter = 'latex')
xticks = 1:3;
xticklabels = ["40", "45", "50"];
ylabel('SNR at gaps [dB]', Interpreter = 'latex')
ax = gca;
ax.FontSize = 15;
ax.TickLabelInterpreter = 'latex';
