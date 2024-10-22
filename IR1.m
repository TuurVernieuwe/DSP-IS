% Channel estimation according to the IR1 method

%% Cleanup
clear; clc; close all

%% Initialize script parameters.
fs = 16000; % Sampling frequency [Hz]
dftsize = 528; % Discrete Fourier Transform (DFT) size [Samples] 
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
% filter = floor(80*abs(h));
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

% %% Filtered white noise Vs. recorded white noise.
% sig = wgn(); % Generate white noise signal    
% 
% % Convolve the IR with the white noise
% gen_noise = fftfilt();                  
% 
% % Send the white noise across the channel
% [] = initparams();
% sim('recplay');
% rec_noise=simout.signals.values(:,1);
% 
% %% Spectrogram and PSD plot
% % Spectrogram
% spectrogram(); % Generated white noise
% spectrogram(); % Recorded white noise
% 
% % PSD
% [] = spectrogram(); % Generated white noise
% [] = spectrogram(); % Recorded white noise
% 
% PSD_gen_noise = ; % Generated white noise PSD
% PSD_rec_noise = ; % Recorded white noise PSD (ideally only when the signal is active)
% 
% figure;
% plot(,pow2db());
% xlabel('Frequency (kHz)');
% ylabel('Power/frequency (dB/Hz)')
% title('Generated white noise.')
% % Output signal
% subplot(2,2)
% plot(,pow2db());
% xlabel('Frequency (kHz)');
% ylabel('Power/frequency (dB/Hz)')
% title('Output signal.')
% sgtitle('Recorded white noise')