% Transmission of a picture across the channel
%% Cleanup
clear; close all; clc;

%% Parameters.
N = 1024; % Total number of symbols in a single OFDM frame, i.e., the DFT size
Lcp = 600; % Cyclic prefix length [samples].
M = 16; %  constellation size.
SNR = 35; % SNR of transmission [dB].
Lt = 20; % Number of training frames.
Ld = 20; % Number of data frames.
fs = 16000; % Sampling frequency [Hz].
channel = 0; % acoustic or simulation

% Bitloading
ON_OFF_mask = ones(N/2-1, 1); % Default all bins to one for regular transmission
bitloading_flag = 0; % If 1 then on-off/adaptive bitloading is enabled.
bitloading_type = "on-off"; % on-off or adaptive 
BW_usage = 1; % Fraction of bins to use for on-off bitloading
Nswitch = (Lt+Ld)*(N+Lcp); % The simulated channel changes every Nswitch number of samples.
smoothing_factor = .99; % Smoothing factor for simulated channel (see simulate_channel.m)

%% Determine bit loading
% Channel estimation based on a dummy transmission for bitloading
train_bits = randi([0 1], log2(M)*(N/2-1), 1); % Generate a random vector of bits
trainblock = qam_mod(train_bits, M); % QAM modulate

if bitloading_flag
    % OFDM modulation of the trainblock
    trainStream = ofdm_mod();

    % Dummy transmission
    if channel == "simulation"
        aligned_Rx = simulate_channel(trainStream, Nswitch,'channel_session6.mat',smoothing_factor);
        aligned_Rx = awgn(aligned_Rx,SNR,'measured');
    else
        [] = initparams();
        size(simin)
        sim('recplay');
        Rx = simout.signals.values(:,1);
        aligned_Rx = alignIO();
    end

    % Extract the estimated channels
    [~, CHANNELS] = ofdm_demod();

    if bitloading_type == "on-off"
        % Only keep bins of high energy
        ON_OFF_mask() = ; % ON-OFF mask with 1 denoting the usage of a bin.
    elseif bitloading_type == "adaptive"
        M = ;     % Constellation sizes
    end
end

%% Construct QAM symbol stream.
% Data blocks
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');
qamStream = qam_mod(bitStream, M);

%% OFDM modulation
[ Tx, nbPackets ] = ofdm_mod( );

%% Transmit OFDM sequence.
if channel == "simulation"
    aligned_Rx = simulate_channel(Tx, Nswitch,'channel_session6.mat',smoothing_factor);
    aligned_Rx = awgn(aligned_Rx,SNR,'measured');
elseif channel == "acoustic"
    [] = initparams();
    size(simin)
    sim('recplay');
    Rx = simout.signals.values(:,1);
    aligned_Rx = alignIO( );
end

%% OFDM Demodulate
[rx_qam, CHANNELS] = ofdm_demod();

%% QAM Demodulate
rx_bits = qam_demod();

%% Bit error rate
BER = ber( )