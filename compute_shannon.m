%% Cleanup
clear; clc; close all

%% Step 1: Initialization
% Hardcoded parameters
fs = 44100;                 % Sampling frequency [Hz]
nyquist_bandwidth = fs/2;   % nyquist bandwidth
N = 1024;                   % DFT size, number of frequency bins
K = N / 2;                  % Half of the DFT bins

% Duration of the recordings
duration = 5;  % Duration for recording signal and noise in seconds

%% Step2: record noise
sig = zeros(duration+fs, 1);
[simin, nbsecs, fs] = initparams(sig, fs);

% Call to recplay.mdl to play simin and record simout
sim('recplay'); 

% Retrieve recorded output
out=simout.signals.values(:,1); 

% get PSD of recorded signal
[pxx, ~] = pwelch(out, hamming(N), K, N, fs);
P_noise = pxx;

%% step3: record signal + noise
sig = wgn(duration * fs, 1, 1);
[simin, nbsecs, fs] = initparams(sig, fs);

% Call to recplay.mdl to play simin and record simout
sim('recplay'); 

% Retrieve recorded output
out=simout.signals.values(:,1); 

% get PSD of recorded signal
[pxx, ~] = pwelch(out, hamming(N), K, N, fs);
P_signalnoise = pxx;

%% Step 4: Compute signal power by subtracting noise PSD
P_signal = P_signalnoise - P_noise;  % Signal PSD is the difference of signal+noise and noise PSDs
P_signal(P_signal < 0) = 0; % Ensure no negative values in PSD due to numerical errors

%% Step 5: Compute the Shannon capacity
Cchannel = 0;   % Initialize channel capacity
for k = 1:K
   Cchannel = Cchannel + log2(1 + P_signal(k) / P_noise(k));  % Sum the contributions from each frequency bin
end
Cchannel = (fs / N) * Cchannel;  % Apply the scaling factor fs / N