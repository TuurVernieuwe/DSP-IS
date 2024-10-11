% Plays and records signals, which are subsequently analyzed by looking at
% the spectrograms, and the power spectral densities.

%% Cleanup
clear; close all;

%% Initialize script parameters
fs = ; % Sampling frequency [Hz]
N = ; % Discrete Fourier Transform (DFT) size [Samples] 
% Overlap length between subsequent frames in the
% short-time-Fourier-transform (STFT) used to plot the spectrogra [samples]
Noverlap = ; 

%% Construct signals
sig = ; 

%% Play and record.
% Call to initparams()
[] = initparams();
% Call to recplay.mdl to play simin and record simout
sim('recplay');
% Retrieve recorded output
out=simout.signals.values(:,1);

%% Compute and plot the spectrogram
% Input signal
figure; subplot(2,1)
spectrogram()
title('Input signal.')

% Output signal
subplot(2,2)
spectrogram()
title('Output signal.')

%% Select input and output signals to compute PSD
% Ideally the output signal should be trimmed such that only periods where
% the signal is active are considered as the PSD assumes stationarity. 
% Else, the PSD will be biased due to the inclusion of the silence.
in = ;
out = ;

%% Compute and plot the power spectral density (PSD)...
% ...Using Welch's method
% Input signal
[] = spectrogram(,,,N,fs);
PSD_Welch_input = ;
% Output signal
[] = spectrogram(,,,N,fs);
PSD_Welch_output = ;
% Plot results
figure;
% Input signal
figure; subplot(2,1)
plot(,pow2db());
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Input signal.')
% Output signal
subplot(2,2)
plot(,pow2db());
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Output signal.')
sgtitle('Power Spectral Density estimate using Welch''s method')

% ...Using Bartlett's method
% Input signal
[] = spectrogram(,rectwin(),,N,fs);
PSD_Bartlett_input = ;
% Output signal
[] = spectrogram(,rectwin(),,N,fs);
PSD_Bartlett_output = ;
% Plot results
figure;
% Input signal
figure; subplot(2,1)
plot(,pow2db());
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Input signal.')
% Output signal
subplot(2,2)
plot(,pow2db());
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Output signal.')
sgtitle('Power Spectral Density estimate using Bartlett''s method')

% ...using the magnitude squared of the frequency spectrum
% Input signal
Power_fft_squared_input = ; % Compute magnitude absolute value squared of the fft of the input
% Rescaling and computations to calcualte the one sided PSD to be
% consistent with 'periodogram' (You can ignore the following lines)
Power_fft_squared_input = Power_fft_squared_input/(length(in)*fs);
PSD_fft_squared_input = Power_fft_squared_input(1:floor(length(in)/2)+1);
PSD_fft_squared_input(2:ceil(length(in)/2)) = ...
    PSD_fft_squared_input(2:ceil(length(in)/2)) + ...
    flipud(Power_fft_squared_input(floor(length(in)/2)+2:length(in)));
% Ouput signal
Power_fft_squared_output = ; % Compute magnitude absolute value squared of the fft of the output
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
figure; subplot(2,1)
plot((0:fs/length():fs/2),pow2db());
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Input signal.')
% Output signal
subplot(2,2)
plot((0:fs/length():fs/2),pow2db());
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Output signal.')
sgtitle('Power Spectral Density estimate using the magnitude squared of the frequency spectrum')