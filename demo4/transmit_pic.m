% Transmission of a picture across the channel
%% Cleanup
clear; close all; clc;

%% Parameters
% Exercise session 7: DMT-OFDM transmission with DD LMS channel tracking.
N = ; % Total number of symbols in a single OFDM frame, i.e., the DFT size
Lcp = ; % Cyclic prefix length [samples].
M = ; % QAM constellation size.
SNR = ; % SNR of transmission [dB].
Lt = ; % Number of training frames.
fs = ; % Sampling frequency [Hz].
channel = ; % simulation or acoustic
mu = ; % NLMS stepsize
alpha = ; % NLMS regularization
Nswitch = (Lt)*(N+Lcp); % The simulated channel changes every Nswitch number of samples.
smoothing_factor = .99; % Smoothing factor for simulated channel (see simulate_channel.m)

% bitloading
ON_OFF_mask = ; % Default all bins to one for regular transmission
BW_usage = ; % Fraction of bins to use for on-off bitloading
bitloading_flag = ; % If 1 then on-off bitloading is enabled.

%% Determine bit loading
% Channel estimation.
% Construct train block.
train_bits = randi(); % Generate a random vector of bits
trainblock = qam_mod(); % QAM modulate 

if bitloading_flag
    % sending only trainblock through the channel to estimate channel    
    trainStream = ofdm_mod( );
    
    if channel == "simulation"
        aligned_Rx = simulate_channel(trainStream, Nswitch,'channel_session7.mat',smoothing_factor); 
        aligned_Rx = awgn(aligned_Rx,SNR,'measured');
    else
        [] = initparams();
        size(simin)
        sim('recplay');
        Rx = simout.signals.values(:,1);
        aligned_Rx = alignIO();
    end
    
    % Compute the channel power    
    [~, CHANNELS] = ofdm_demod( );

    % Only keep bins of high energy
    ON_OFF_mask() = ;
end

%% Construct QAM symbol stream.
% Data blocks
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');
qamStream = qam_mod();

%% OFDM modulation.
Tx = ofdm_mod();

%% Transmit OFDM sequence.
if channel == "simulation"
    aligned_Rx = simulate_channel(Tx, Nswitch,'channel_session7.mat',smoothing_factor);
    aligned_Rx = awgn(aligned_Rx,SNR,'measured');
elseif channel == "acoustic"
    [] = initparams();
    sim('recplay');
    Rx = simout.signals.values(:,1);
    aligned_Rx = alignIO( );
end

%% OFDM Demodulate
[rx_qam, CHANNELS] = ofdm_demod();

%% QAM demodulation.
rx_bits = qam_demod();

%% Bit error rate
BER = ber( )