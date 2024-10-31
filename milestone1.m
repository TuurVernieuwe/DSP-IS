%% Cleanup
clear; close all;

%%zeek1
% Initialize script parameters
fs = 16000; % Sampling frequency [Hz]
N = 512; % Discrete Fourier Transform (DFT) size [Samples] 
t_max = 2; % Length of the signal [s]
t = 0:1/fs:t_max;
f0 = 1500;
% Overlap length between subsequent frames in the
% short-time-Fourier-transform (STFT) used to plot the spectrogra [samples]
Noverlap = 256; 

% Construct signals
sig = wgn(2*fs, 1, 1);

% Play and record.
% Call to initparams()
[simin, nbsecs, fs] = initparams(sig, fs);
% Call to recplay.mdl to play simin and record simout
sim('recplay');
% Retrieve recorded output
out=simout.signals.values(:,1);


% Compute and plot the spectrogram
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
title('Spectrogram of Signal out "out"');
colorbar;

% Select input and output signals to compute PSD
% Ideally the output signal should be trimmed such that only periods where
% the signal is active are considered as the PSD assumes stationarity. 
% Else, the PSD will be biased due to the inclusion of the silence.
in = simin(2*fs:end-fs, 1);
out_selected = out(2*fs:end-fs, 1);

% Compute and plot the power spectral density (PSD)...
% ...Using Welch's method
% Input signal
figure; subplot(2,1,1);
[pxx_sig, f_sig] = pwelch(in, N, Noverlap, N, fs);  % Compute PSD using Welch's method
plot(f_sig, 10*log10(pxx_sig));  % Convert to dB scale
title('PSD of Transmitted Signal (Welch Method)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
ylim([-54, -48]);
grid on;
% Output signal
subplot(2,1,2);
[pxx_mic, f_mic] = pwelch(out_selected, N, Noverlap, N, fs);  % Compute PSD using Welch's method
plot(f_mic, 10*log10(pxx_mic));  % Convert to dB scale
title('PSD of Recorded Signal (Welch Method)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
ylim([-100, -50]);
grid on;

% ...Using Bartlett's method
% Input signal
[S,F,T,P] = spectrogram(in,N,0,N,fs);
PSD_Bartlett_input = mean(P, 2);
% Output signal
[S,F_mic,T,P] = spectrogram(out,N,0,N,fs);
PSD_Bartlett_output = mean(P, 2);
% Plot results
% Input signal
figure; subplot(2,1,1)
plot(F,pow2db(PSD_Bartlett_input));
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('PSD of Transmitted Signal (Barlett Method)');
ylim([-54, -48]);
grid on;
% Output signal
subplot(2,1,2)
plot(F_mic,pow2db(PSD_Bartlett_output));
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('PSD of Recorded Signal (Barlett Method)');
ylim([-100, -50]);
grid on;

%% week2: IR2
% Channel estimation according to the IR2 method

% Initialize script parameters.
dftsize = 256; % Discrete Fourier Transform (DFT) size [Samples] 
% Overlap length between subsequent frames in the
% short-time-Fourier-transform (STFT) [samples]
channelLength = 550; % Length of impulse response [samples]
delay = 6720; % Positive delay safety margin when aligning input and output [samples]

% Calculate the impulse response.
uMatrix = toeplitz(sig , [sig(1); zeros(channelLength-1, 1)]); % Toeplitz matrix
yOnset = 2*fs+delay; % Determine start of recorded signal [samples]
y = out(yOnset:yOnset+size(uMatrix, 1)-1); % Extract the relevant output signal

h = lsqr(uMatrix, y); % Estimate impulse response

save('channel.mat','h'); % Save impulser response

% Plot IR.
% Time domain signal
figure; subplot(2,1,1)
plot(h);
xlabel('Time [samples]')
ylabel('Impulse response [arb.]')
title("Time response")
% Magnitude response
subplot(2,1,2)
plot(pow2db(abs(fft(h, dftsize))));
xlabel('Frequency [Hz]')
ylabel('Magnitude response [dB]')
title("Frequency response")
sgtitle("IR2")

%% week2: IR1
%% Create the signal to be played
sig = [1; zeros(2*fs-1, 1)];                       

%% Play and record.
% Call to initparams()
[simin, nbsecs, fs] = initparams(sig, fs);
% Call to recplay.mdl to play simin and record simout
sim('recplay');
% Retrieve recorded output
out=simout.signals.values(:,1);

%% Calculate the IR by trimming the output
h = out(2*fs:end-fs, 1);
filter = floor(1000*abs(h));
h = filter.*h;

%% Plot IR.
% Time domain signal
figure; subplot(2,1,1)
plot(h);
xlabel('Time [samples]')
ylabel('Impulse response [arb.]')
title("Time response")
% Magnitude response
subplot(2,1,2)
plot(pow2db(abs(fft(h, dftsize))));
xlabel('Frequency [Hz]')
ylabel('Magnitude response [dB]')
title("Frequency response")
sgtitle("IR1")
