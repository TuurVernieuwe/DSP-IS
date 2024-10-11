% Plays and records signals, which are subsequently analyzed by looking at
% the spectrograms, and the power spectral densities.

%% Cleanup
clear; close all;

%% Initialize script parameters
fs = 16000; % Sampling frequency [Hz]
t_max = 2; % Length of the signal [s]
N = 10; % Discrete Fourier Transform (DFT) size [Samples] 
% Overlap length between subsequent frames in the
% short-time-Fourier-transform (STFT) used to plot the spectrogra [samples]
Noverlap = 5; 

%% Construct signals
sig = sin(2*pi*linspace(0, 400*t_max, fs*t_max)'); 

%% Play and record.
% Call to initparams()
[simin, nbsecs, fs] = initparams(sig, fs);
% Call to recplay.mdl to play simin and record simout
sim('recplay');
% Retrieve recorded output
out=simout.signals.values(:,1);

%% Compute and plot the spectrogram
% Input signal
figure; subplot(2,1,1)
spectrogram(sig, N, Noverlap)
title('Input signal.')
xlabel('Time (s)'); ylabel('Frequency (Hz)')

% Output signal
subplot(2,1,2)
spectrogram(out, N, Noverlap)
title('Output signal.')
xlabel('Time (s)'); ylabel('Frequency (Hz)')

%% Select input and output signals to compute PSD
% Ideally the output signal should be trimmed such that only periods where
% the signal is active are considered as the PSD assumes stationarity. 
% Else, the PSD will be biased due to the inclusion of the silence.
in = simin(2*fs:end-fs, 1);
out = out(2*fs:end-fs, 1);

%% Compute and plot the power spectral density (PSD)...
% ...Using Welch's method
% Input signal
[S,F,T,P] = spectrogram(in,N,Noverlap,N,fs);
PSD_Welch_input = P;
% Output signal
[S,F,T,P] = spectrogram(out,N,Noverlap,N,fs);
PSD_Welch_output = P;
% Plot results
figure;
% Input signal
figure; subplot(2,1,1)
plot(F,pow2db(PSD_Welch_input));
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Input signal.')
% Output signal
subplot(2,1,2)
plot(F,pow2db(PSD_Welch_output));
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Output signal.')
sgtitle('Power Spectral Density estimate using Welch''s method')

% ...Using Bartlett's method
% Input signal
[S,F,T,P] = spectrogram(sig,rectwin(),Noverlap,N,fs);
PSD_Bartlett_input = P;
% Output signal
[S,F,T,P] = spectrogram(out,rectwin(),Noverlap,N,fs);
PSD_Bartlett_output = P;
% Plot results
figure;
% Input signal
figure; subplot(2,1,1)
plot(F,pow2db(PSD_Bartlett_input));
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Input signal.')
% Output signal
subplot(2,1,2)
plot(F,pow2db(PSD_Bartlett_output));
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Output signal.')
sgtitle('Power Spectral Density estimate using Bartlett''s method')

% ...using the magnitude squared of the frequency spectrum
% Input signal
Power_fft_squared_input = abs(fft(sig))^2; % Compute magnitude absolute value squared of the fft of the input
% Rescaling and computations to calcualte the one sided PSD to be
% consistent with 'periodogram' (You can ignore the following lines)
Power_fft_squared_input = Power_fft_squared_input/(length(in)*fs);
PSD_fft_squared_input = Power_fft_squared_input(1:floor(length(in)/2)+1);
PSD_fft_squared_input(2:ceil(length(in)/2)) = ...
    PSD_fft_squared_input(2:ceil(length(in)/2)) + ...
    flipud(Power_fft_squared_input(floor(length(in)/2)+2:length(in)));
% Ouput signal
Power_fft_squared_output = abs(fft(out))^2; % Compute magnitude absolute value squared of the fft of the output
% Rescaling and computations to calculate the one sided PSD to be
% consistent with 'periodogram' (You can ignore the following lines)
Power_fft_squared_output = Power_fft_squared_output/(length(out)*fs);
PSD_fft_squared_output = Power_fft_squared_output(1:floor(length(out)/2)+1);
PSD_fft_squared_output(2:ceil(length(out)/2)) = ...
    PSD_fft_squared_output(2:ceil(length(out)/2)) + ...
    flipud(Power_fft_squared_output(floor(length(out)/2)+2:length(out)));

% Plot results
figure;
% Input signal
figure; subplot(2,1,1)
plot((0:fs/length(PSD_fft_squared_input):fs/2),pow2db(PSD_fft_squared_input));
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Input signal.')
% Output signal
subplot(2,1,2)
plot((0:fs/length(PSD_fft_squared_output):fs/2),pow2db(PSD_fft_squared_output));
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Output signal.')
sgtitle('Power Spectral Density estimate using the magnitude squared of the frequency spectrum')