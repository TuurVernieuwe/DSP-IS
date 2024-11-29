% Transmission of a picture across the channel
%% Cleanup
clear; close all; clc;

%% Parameters.
N = 1024; % Total number of symbols in a single OFDM frame, i.e., the DFT size
Lcp = fix(N/4); % Cyclic prefix length [samples].
M = 16; %  constellation size.
SNR = 15; % SNR of transmission [dB].
Lt = 5; % Number of training frames.
Ld = 5; % Number of data frames.
fs = 16000; % Sampling frequency [Hz].
channel = "simulation"; % acoustic or simulation

% Bitloading
ON_OFF_mask = ones(N/2-1, 1); % Default all bins to one for regular transmission
bitloading_flag = 1; % If 1 then on-off/adaptive bitloading is enabled.
bitloading_type = "on-off"; % on-off or adaptive 
BW_usage = 0.6; % Fraction of bins to use for on-off bitloading
Nswitch = (Lt+Ld)*(N+Lcp); % The simulated channel changes every Nswitch number of samples.
smoothing_factor = .99; % Smoothing factor for simulated channel (see simulate_channel.m)

% Acoustic
f_sync = 5000;
t = 0:1/fs:1;

%% Determine bit loading
% Channel estimation based on a dummy transmission for bitloading
train_bits = randi([0 1], log2(M)*sum(ON_OFF_mask), 1); % Generate a random vector of bits
trainblock = qam_mod(train_bits, M); % QAM modulate

if bitloading_flag
    % OFDM modulation of the trainblock
    trainStream = ofdm_mod(repmat(trainblock, 10, 1), N, Lcp, ON_OFF_mask);

    % Dummy transmission
    if channel == "simulation"
        aligned_Rx = simulate_channel(trainStream, Nswitch,'channel_session6.mat',smoothing_factor);
        aligned_Rx = awgn(aligned_Rx,SNR,'measured');
    else
        [simin, nbsecs, fs] = initparams(trainStream, fs, sin(2*pi*t*f_sync)');
        sim('recplay');
        Rx = simout.signals.values(:,1);
        aligned_Rx = alignIO(Rx, fs);
        aligned_Rx = aligned_Rx(1:length(trainStream));
    end

    % Extract the estimated channels
    [~, CHANNELS] = ofdm_demod(aligned_Rx, N, Lcp, 10*length(trainblock), ON_OFF_mask, trainblock);

    if bitloading_type == "on-off"
        % Only keep bins of high energy
        sorted = sort(abs(CHANNELS), 'descend');
        idx = max(floor(BW_usage*length(CHANNELS)), 1);
        threshold = sorted(idx);
        ON_OFF_mask = abs(CHANNELS) >= threshold; % ON-OFF mask with 1 denoting the usage of a bin.
        trainblock = trainblock(logical(ON_OFF_mask));
    elseif bitloading_type == "adaptive"
        M = 0; % Constellation sizes
    end
end

%% Construct QAM symbol stream.
% Data blocks
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');
qamStream = qam_mod(bitStream, M);

%% OFDM modulation
[ Tx, nbPackets ] = ofdm_mod(qamStream, N, Lcp, ON_OFF_mask, Lt, Ld, trainblock);

%% Transmit OFDM sequence.
if channel == "simulation"
    aligned_Rx = simulate_channel(Tx, Nswitch,'channel_session6.mat',smoothing_factor);
    aligned_Rx = awgn(aligned_Rx,SNR,'measured');
elseif channel == "acoustic"
    [simin,nbsecs,fs] = initparams(Tx, fs, sin(2*pi*t*f_sync)');
    sim('recplay');
    Rx = simout.signals.values(:,1);
    aligned_Rx = alignIO(Rx, fs);
end

%% OFDM Demodulate
tic
[rx_qam, CHANNELS] = ofdm_demod(aligned_Rx, N, Lcp, length(qamStream), ON_OFF_mask, trainblock, Lt, Ld, nbPackets);
time = toc;
time_estimate = time/nbPackets;

%% QAM Demodulate
rx_bits = qam_demod(rx_qam, M, length(bitStream));

%% Bit error rate
BER = ber(bitStream, rx_bits)