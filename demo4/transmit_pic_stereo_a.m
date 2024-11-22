% Transmission of a picture across the channel
%% Cleanup
clear; close all; clc;

% Parameters.
Lh = ; % Length of impulse response
N = ; % Total number of symbols in a single OFDM frame, i.e., the DFT size
M = ; % QAM constellation size.
Lcp = ; % Cyclic prefix length [samples].
Lt = ; % Number of training frames
Equalization = "adaptive"; % Equalization mode (see ofdm_demod_stereo.m)
mu = ;% NLMS stepsize
alpha = ; % NLMS regularization
SNR = ; % SNR of transmission 

%% Generate two random impulse responses, and calculate frequency response.
h1 = ; h2 = ; % Impulse responses
H = ; % N/2-1X2 matrix containing frequency transform of h1 and h2

%% Construct QAM symbol stream.
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');
qamStream = qam_mod( );

% Construct train block.
train_bits = randi();
train_frame = qam_mod();

%% OFDM modulation.
[Tx, a, b] = ofdm_mod_stereo();

%% Transmit symbol.
Rx = ;
Rx = awgn();

%% OFDM demodulation.
[rec_qamStream, CHANNELS] = ofdm_demod_stereo();

%% QAM demodulation.
rx_bits = qam_demod();

%% Calculate BER
BER = ber( )