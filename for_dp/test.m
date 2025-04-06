% % Parameters
% N = 16; % Window length
% window = hann(N); % Generate the Hann window
% 
% % Compute the Fourier Transform
% fft_window = fft(window, 2048); % Zero-padding for better frequency resolution
% magnitude_spectrum = 20*log(abs(fftshift(fft_window))); % Center the spectrum
% 
% % Frequency axis
% fs = 1; % Assume a normalized sampling rate
% f = linspace(-fs/2, fs/2, length(magnitude_spectrum));
% 
% % Show only specific portions of the spectrum
% portion_to_show = (f > -20 & f < 20); % Change this range as needed
% 
% % Plotting
% figure;
% plot(f, magnitude_spectrum, 'b', 'LineWidth', 1.5);
% grid on;
% xlabel('Frequency (normalized)');
% ylabel('Magnitude');
% title('Fourier Transform of a Hann Window');
% Trim or pad to the target number of samples
% Load the input .wav file
[inputFile, fs] = audioread('II. Double.mp3'); % Replace 'input.wav' with your file's name

% Set target parameters
targetDuration = 4; % seconds
targetSampleRate = 16000; % Hz
targetSamples = targetSampleRate * targetDuration;

% Resample the audio if necessary
if fs ~= targetSampleRate
    inputFile = resample(inputFile, targetSampleRate, fs);
    fs = targetSampleRate; % Update the sample rate
end

% Convert to mono if stereo
if size(inputFile, 2) == 2
    inputFile = mean(inputFile, 2); % Take the average of the two channels
end

% Trim or pad to the target number of samples
currentSamples = length(inputFile);
if currentSamples > targetSamples
    % Trim the audio
    outputAudio = inputFile(1:targetSamples);
elseif currentSamples < targetSamples
    % Pad the audio with zeros
    outputAudio = [inputFile; zeros(targetSamples - currentSamples, 1)];
else
    % No change needed
    outputAudio = inputFile;
end

% Save the processed audio to a new file
audiowrite('violin_input.wav', outputAudio, targetSampleRate);

disp('Audio has been processed and saved as "violioninput.wav".');
