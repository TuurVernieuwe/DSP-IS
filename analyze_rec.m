% Plays and records signals, which are subsequently analyzed by looking at
% the spectrograms, and the power spectral densities.

%% Cleanup
clear; close all;

%% Initialize script parameters
fs = 16000; % Sampling frequency [Hz]
N = 256; % Discrete Fourier Transform (DFT) size [Samples] 
t_max = 2; % Length of the signal [s]
t = 0:1/fs:t_max;
f0 = 400;
% Overlap length between subsequent frames in the
% short-time-Fourier-transform (STFT) used to plot the spectrogra [samples]
Noverlap = 128; 

%% Construct signals
sig = wgn(2*fs,1,0); 

%% Play and record.
% Call to initparams()
[simin, nbsecs, fs] = initparams(sig, fs);
% Call to recplay.mdl to play simin and record simout
sim('recplay');
% Retrieve recorded output
out=simout.signals.values(:,1);

%% Compute and plot the spectrogram
% Input signal
figure; subplot(2,1,1);        % First subplot for the signal in 'sig'
[S, F, T, P] = spectrogram(sig, N, Noverlap, N, fs);   % Compute spectrogram
imagesc(T, F, 10*log10(P));                            % Plot in dB scale
axis xy;                                                % Correct axis orientation
xlabel('Time (s)');
ylabel('Frequency (Hz)');
ylim([0 2000]);
title('Spectrogram of Signal in "sig"');
colorbar;

% Output signal
subplot(2,1,2);        % First subplot for the signal in 'sig'
[S_mic, F_mic, T_mic, P_mic] = spectrogram(out, N, Noverlap, N, fs);   % Compute spectrogram
imagesc(T_mic, F_mic, 10*log10(P_mic));                            % Plot in dB scale
axis xy;                                                % Correct axis orientation
xlabel('Time (s)');
ylabel('Frequency (Hz)');
ylim([0 2000]);
xlim([2, 4])
title('Spectrogram of Signal out "out"');
colorbar;

%% Select input and output signals to compute PSD
% Ideally the output signal should be trimmed such that only periods where
% the signal is active are considered as the PSD assumes stationarity. 
% Else, the PSD will be biased due to the inclusion of the silence.
in = simin(2*fs:end-fs, 1);
out = out(2*fs:end-fs, 1);

%% Compute and plot the power spectral density (PSD)...
% ...Using Welch's method
% Input signal
figure; subplot(2,1,1);
[pxx_sig, f_sig] = pwelch(in, N, Noverlap, N, fs);  % Compute PSD using Welch's method
plot(f_sig, 10*log10(pxx_sig));  % Convert to dB scale
title('PSD of Transmitted Signal (Welch Method)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
grid on;
% Output signal
subplot(2,1,2);
[pxx_mic, f_mic] = pwelch(out, N, Noverlap, N, fs);  % Compute PSD using Welch's method
plot(f_mic, 10*log10(pxx_mic));  % Convert to dB scale
title('PSD of Recorded Signal (Welch Method)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
grid on;

% ...Using Bartlett's method
% Input signal
[S,F,T,P] = spectrogram(in,N,Noverlap,N,fs);
PSD_Bartlett_input = mean(P, 2);
% Output signal
[S,F_mic,T,P] = spectrogram(out,N,Noverlap,N,fs);
PSD_Bartlett_output = mean(P, 2);
% Plot results
% Input signal
figure; subplot(2,1,1)
plot(F,pow2db(PSD_Bartlett_input));
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Input signal.')
grid on;
% Output signal
subplot(2,1,2)
plot(F_mic,pow2db(PSD_Bartlett_output));
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Output signal.')
sgtitle('Power Spectral Density estimate using Bartlett''s method')
grid on;

% ...using the magnitude squared of the frequency spectrum
% Input signal
Power_fft_squared_input = abs(fft(in)).^2; % Compute magnitude absolute value squared of the fft of the input
% Rescaling and computations to calcualte the one sided PSD to be
% consistent with 'periodogram' (You can ignore the following lines)
Power_fft_squared_input = Power_fft_squared_input/(length(in)*fs);
PSD_fft_squared_input = Power_fft_squared_input(1:floor(length(in)/2)+1);
PSD_fft_squared_input(2:ceil(length(in)/2)) = ...
    PSD_fft_squared_input(2:ceil(length(in)/2)) + ...
    flipud(Power_fft_squared_input(floor(length(in)/2)+2:length(in)));
% Ouput signal
Power_fft_squared_output = abs(fft(out)).^2; % Compute magnitude absolute value squared of the fft of the output
% Rescaling and computations to calculate the one sided PSD to be
% consistent with 'periodogram' (You can ignore the following lines)
Power_fft_squared_output = Power_fft_squared_output/(length(out)*fs);
PSD_fft_squared_output = Power_fft_squared_output(1:floor(length(out)/2)+1);
PSD_fft_squared_output(2:ceil(length(out)/2)) = ...
    PSD_fft_squared_output(2:ceil(length(out)/2)) + ...
    flipud(Power_fft_squared_output(floor(length(out)/2)+2:length(out)));

% Plot results
% Input signal
figure; subplot(2,1,1)
plot((0:fs/length(Power_fft_squared_input):fs/2),pow2db(PSD_fft_squared_input));
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Input signal.')
% Output signal
subplot(2,1,2)
plot((0:fs/length(Power_fft_squared_output):fs/2),pow2db(PSD_fft_squared_output));
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Output signal.')
sgtitle('Power Spectral Density estimate using the magnitude squared of the frequency spectrum')