% Transmission of a picture across the channel
%% Cleanup
clear; close all; clc;

% Parameters.
Lh = 200; % Length of impulse response
N = 1024; % Total number of symbols in a single OFDM frame, i.e., the DFT size
M = 16; % QAM constellation size.
Lcp = 300; % Cyclic prefix length [samples].
Lt = 10; % Number of training frames
Ld = 5; % Does not matter in this file
Equalization = "adaptive"; % Equalization mode (see ofdm_demod_stereo.m)
mu = 0.2;% NLMS stepsize
alpha = 1; % NLMS regularization
SNR = inf; % SNR of transmission 

%% Generate two random impulse responses, and calculate frequency response.
h1 = randn(Lh, 1); h2 = randn(Lh, 1); % Impulse responses
H = [fft(h1, N), fft(h2, N)]; % N/2-1X2 matrix containing frequency transform of h1 and h2

%% Construct QAM symbol stream.
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');
qamStream = qam_mod(bitStream, M);

% Construct train block.
trainbits = randi([0, 1], log2(M)*(N/2-1), 1);
trainblock = qam_mod(trainbits, M);

%% OFDM modulation.
[Tx, a, b, nbpackets] = ofdm_mod_stereo(qamStream, N, Lcp, H, Lt, Ld, trainblock, 'fixed');

%% Transmit symbol.
Rx = fftfilt(h1, Tx(:,1)) + fftfilt(h2, Tx(:,2));
Rx = awgn(Rx, SNR, 'measured');

%% OFDM demodulation.
[rec_qamStream, CHANNELS] = ofdm_demod_stereo(Rx, N, Lcp, trainblock, Lt, Ld, M, nbpackets, 'fixed', mu, alpha);

%% QAM demodulation.
rx_bits = qam_demod(rec_qamStream, M, length(bitStream));

%% Calculate BER
BER = ber(bitStream, rx_bits)