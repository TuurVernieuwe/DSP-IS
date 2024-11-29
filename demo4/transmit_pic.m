% Transmission of a picture across the channel
%% Cleanup
clear; close all; clc;

%% Parameters
% Exercise session 7: DMT-OFDM transmission with DD LMS channel tracking.
N = 1024; % Total number of symbols in a single OFDM frame, i.e., the DFT size
Lcp = fix(N/4); % Cyclic prefix length [samples].
M = 16; % QAM constellation size.
SNR = 30; % SNR of transmission [dB].
Lt = 5; % Number of training frames.
fs = 16000; % Sampling frequency [Hz].
channel = 'simulation'; % simulation or acoustic
mu = 0.2; % NLMS stepsize
alpha = 1; % NLMS regularization
Nswitch = (Lt)*(N+Lcp); % The simulated channel changes every Nswitch number of samples.
smoothing_factor = .99; % Smoothing factor for simulated channel (see simulate_channel.m)

% bitloading
ON_OFF_mask = ones(N/2-1, 1); % Default all bins to one for regular transmission
BW_usage = 1; % Fraction of bins to use for on-off bitloading
bitloading_flag = 0; % If 1 then on-off bitloading is enabled.

% Acoustic
f_sync = 5000;
t = 0:1/fs:1;

%% Determine bit loading
% Channel estimation.
% Construct train block.
train_bits = randi([0 1], log2(M)*sum(ON_OFF_mask), 1); % Generate a random vector of bits
trainblock = qam_mod(train_bits, M); % QAM modulate 

if bitloading_flag
    % sending only trainblock through the channel to estimate channel    
    trainStream = ofdm_mod(trainblock, N, Lcp, ON_OFF_mask);

    if channel == "simulation"
        aligned_Rx = simulate_channel(trainStream, Nswitch,'channel_session7.mat',smoothing_factor); 
        aligned_Rx = awgn(aligned_Rx,SNR,'measured');
    else
        [simin, nbsecs, fs] = initparams(trainStream, fs, sin(2*pi*t*f_sync)');
        size(simin)
        sim('recplay');
        Rx = simout.signals.values(:,1);
        aligned_Rx = alignIO(Rx, fs);
        aligned_Rx = aligned_Rx(1:length(trainStream));
    end

    % Compute the channel power    
    [~, CHANNELS] = ofdm_demod(aligned_Rx, N, Lcp, length(trainblock), ON_OFF_mask, trainblock);

    % Only keep bins of high energy
    sorted = sort(abs(CHANNELS), 'descend');
    idx = max(floor(BW_usage*length(CHANNELS)), 1);
    threshold = sorted(idx);
    ON_OFF_mask = abs(CHANNELS) >= threshold; % ON-OFF mask with 1 denoting the usage of a bin.
    trainblock = trainblock(logical(ON_OFF_mask));
end

%% Construct QAM symbol stream.
% Data blocks
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');
qamStream = qam_mod(bitStream, M);

%% OFDM modulation.
Tx = ofdm_mod(qamStream, N, Lcp, Lt, trainblock, ON_OFF_mask);

%% Transmit OFDM sequence.
if channel == "simulation"
    aligned_Rx = simulate_channel(Tx, Nswitch,'channel_session7.mat',smoothing_factor);
    aligned_Rx = awgn(aligned_Rx,SNR,'measured');
elseif channel == "acoustic"
    [simin,nbsecs,fs] = initparams(Tx, fs, sin(2*pi*t*f_sync)');
    sim('recplay');
    Rx = simout.signals.values(:,1);
    aligned_Rx = alignIO(Rx, fs);
end

%% OFDM Demodulate
[rx_qam, CHANNELS] = ofdm_demod(aligned_Rx, N, Lcp, Lt, M, trainblock, ON_OFF_mask, mu, alpha, 'NLMS');

%% QAM demodulation.
rx_bits = qam_demod(rx_qam, M, length(bitStream));

%% Bit error rate
BER = ber(bitStream, rx_bits)