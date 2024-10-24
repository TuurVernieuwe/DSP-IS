% Channel estimation according to the IR2 method

%% Cleanup
clear; clc; close all

% Initialize script parameters.
fs = 16000; % Sampling frequency [Hz]
dftsize = 256; % Discrete Fourier Transform (DFT) size [Samples] 
% Overlap length between subsequent frames in the
% short-time-Fourier-transform (STFT) [samples]
Noverlap = 128; 
channelLength = 550; % Length of impulse response [samples]
delay = 6720; % Positive delay safety margin when aligning input and output [samples]


%% Create the signal to be played.
duration = 2; % Duration of the signal in [s]
sig = wgn(duration*fs, 1, 1);  

%% Filter signal -> only to be used for exercise 2.3 (To this end, also copy the content
%% of the current file to a new IR_bandstop.m file)
filt = fir1(100, [700 3000] / (fs / 2), 'stop');
sig = fftfilt(filt, sig);

%% Play and record.
% Call to initparams()
[simin, nbsecs, fs] = initparams(sig, fs);
% Call to recplay.mdl to play simin and record simout
sim('recplay');
% Retrieve recorded output
out=simout.signals.values(:,1);

%% Calculate the impulse response.
uMatrix = toeplitz(sig , [sig(1); zeros(channelLength-1, 1)]); % Toeplitz matrix
yOnset = 2*fs + delay; % Determine start of recorded signal [samples]
y = out(yOnset:yOnset+size(uMatrix, 1)-1); % Extract the relevant output signal

h = lsqr(uMatrix, y); % Estimate impulse response

save('channel.mat','h'); % Save impulser response

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
sig = wgn(duration*fs, 1, 1); % Generate white noise signal    

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
spectrogram(gen_noise); % Generated white noise
subplot(2, 1, 2);
spectrogram(rec_noise); % Recorded white noise

% PSD
[S_gen, F_gen, T_gen, P_gen] = spectrogram(gen_noise); % Generated white noise
[S_rec, F_rec, T_rec, P_rec] = spectrogram(rec_noise); % Recorded white noise

PSD_gen_noise = mean(P_gen, 2); % Generated white noise PSD
PSD_rec_noise = mean(P_rec, 2); % Recorded white noise PSD (ideally only when the signal is active)

figure;
subplot(2, 1, 1);
plot(F_gen, pow2db(PSD_gen_noise));
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Generated white noise.')
% Output signal
subplot(2, 1, 2);
plot(F_rec ,pow2db(PSD_rec_noise));
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Output signal.')
sgtitle('Recorded white noise')