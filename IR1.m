% Channel estimation according to the IR1 method

%% Cleanup
clear; clc; close all

%% Initialize script parameters.
fs = 12000; % Sampling frequency [Hz]
dftsize = 528; % Discrete Fourier Transform (DFT) size [Samples] 
N = 512;
% Overlap length between subsequent frames in the
% short-time-Fourier-transform (STFT) [samples]
Noverlap = 256; 

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
% filter = floor(60*abs(h));
% h = filter.*h;

%% Plot IR.
% Time domain signal
figure; subplot(2,1,1)
plot(h);
xlabel('Time [samples]')
ylabel('Impulse response [arb.]')
% Magnitude response
subplot(2,1,2)
plot(pow2db(abs(fft(h, dftsize))));
xlabel('Frequency [Hz]')
ylabel('Magnitude response [dB]')

%% Filtered white noise Vs. recorded white noise.
sig = wgn(2*fs, 1, 1); % Generate white noise signal    

% Convolve the IR with the white noise
gen_noise = fftfilt(h, sig);                  

% Send the white noise across the channel
[simin, nbsecs, fs] = initparams(sig, fs);
sim('recplay');
rec_noise=simout.signals.values(:,1);

%% Spectrogram and PSD plot
% Spectrogram
figure;
subplot(2, 1, 1);
spectrogram(gen_noise, N, Noverlap, N, fs); % Generated white noise
subplot(2, 1, 2);
spectrogram(rec_noise, N, Noverlap, N, fs); % Recorded white noise

% PSD
[pxx_gen, f_gen] = pwelch(gen_noise, N, Noverlap, N, fs); % Generated white noise
[pxx_rec, f_rec] = pwelch(rec_noise, N, Noverlap, N, fs); % Recorded white noise

PSD_gen_noise = pxx_gen; % Generated white noise PSD
PSD_rec_noise = pxx_rec; % Recorded white noise PSD (ideally only when the signal is active)

figure;
subplot(2, 1, 1)
plot(f_gen, pow2db(PSD_gen_noise));
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Generated white noise.')
% Output signal
subplot(2, 1, 2)
plot(f_rec, pow2db(PSD_rec_noise));
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Output signal.')
sgtitle('Recorded white noise')