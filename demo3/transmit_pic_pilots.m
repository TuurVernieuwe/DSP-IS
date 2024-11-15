%% Cleanup
clear; clc; close all;

%% Parameters.
N = ; % Total number of symbols in a single OFDM frame, i.e., the DFT size
Lcp = ; % Cyclic prefix length [samples].
M = ; % QAM constellation size.
SNR = ; % SNR of transmission [dB].
fs = ; % Sampling frequency [Hz].
channel = ; % simulation or acoustic
Nswitch = N+Lcp; % The simulated channel changes every Nswitch number of samples.
smoothing_factor = 0.99; % Smoothing factor for simulated channel (see simulate_channel.m)

%% Construct QAM symbol stream.
% Data blocks
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');
qamStream = qam_mod();

%% OFDM modulation
[Tx, nbOFDMsymb] = ofdm_mod_pilots();

%% Transmit OFDM sequence.
if channel == "simulation"
    aligned_Rx = simulate_channel(Tx, Nswitch,'channel_session6.mat',smoothing_factor);
    aligned_Rx = awgn(aligned_Rx,SNR,'measured');
elseif channel == "acoustic"
    [] = initparams();
    size(simin)
    sim('recplay');
    Rx = simout.signals.values(:,1);
    aligned_Rx = alignIO();
end

%% OFDM Demodulate
[rx_qam, CHANNELS] = ofdm_demod_pilots();

%% QAM Demodulate
rx_bits = qam_demod();

%% Bit error rate
BER = ber( )