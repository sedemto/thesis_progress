clc;clear all;
% Parameters
fs = 100;               % Sampling frequency (Hz)
t = linspace(0, 1, fs);       % Time vector (1 second)
signal = randn(size(t)); % Random audio signal
amplitude = 1;              % Random amplitude between 0 and 2
frequency = 10;             % Random frequency up to 10 Hz
frequency2 = 12;
phase = pi/2;             % Random phase between 0 and 2*pi

% Generate the sine wave
sine_wave = amplitude * (sin(2 * pi * frequency * t + phase));
signal = sine_wave;
% Compute STFT
window_length = 8;            % Length of the STFT window
overlap_length = window_length / 2;
nfft = window_length;
[S, f, n] = stft(signal, fs, 'Window', hamming(window_length), ...
                 'OverlapLength', overlap_length, 'FFTLength', nfft);

% Separate magnitude and phase
magnitude = abs(S);
phase = angle(S);
% Scenario 1: Corrupted Phase (phase replaced by zeros)


corrupted_phase = phase(:,:);
corrupted_phase(:,7:18) = 1;

corrupted_phase_signal = istft(magnitude .* exp(1j * corrupted_phase), fs, ...
                               'Window', hamming(window_length), ...
                               'OverlapLength', overlap_length, 'FFTLength', nfft);

% Scenario 2: Corrupted Magnitude (magnitude replaced by ones)
corrupted_magnitude = magnitude(:,:);
corrupted_magnitude(:,7:18) = 1;
corrupted_magnitude_signal = istft(corrupted_magnitude .* exp(1j * phase), fs, ...
                                   'Window', hamming(window_length), ...
                                   'OverlapLength', overlap_length, 'FFTLength', nfft);

% Plotting the original, phase-corrupted, and magnitude-corrupted signals
figure;
plot(t, signal, 'k', 'DisplayName', 'Original Signal'); hold on;
plot(t, real(corrupted_phase_signal), 'r', 'DisplayName', 'Phase Corrupted Signal');
plot(t, real(corrupted_magnitude_signal), 'b', 'DisplayName', 'Magnitude Corrupted Signal');
hold off;

% Add titles, labels, and legend
xlabel('Time (s)');
ylabel('Amplitude');
legend('show');
