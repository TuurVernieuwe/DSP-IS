% Transmission of a picture across the channel
%% Cleanup
clear; close all; clc;

% Parameters.
Lh = 200; % Length of impulse response
N = 1024; % Total number of symbols in a single OFDM frame, i.e., the DFT size
M = 16; % QAM constellation size.
Lcp = Lh; % Cyclic prefix length [samples].
Lt = 5; % Number of training frames
Ld = 5; % Number of data frames
Equalization = "fixed"; % Equalization mode (see ofdm_demod_stereo.m)
mu = 0.2;% NLMS stepsize
alpha = 1; % NLMS regularization
SNR = 20; % SNR of transmission 

%% Generate two random impulse responses, and calculate frequency response.
load('channel_stereo_session7.mat') % Impulse responses
H = [fft(h1(1:min(length(h1),Lh)), N) fft(h2(1:min(length(h1),Lh)), N)];
H = H(1:N/2-1,:); % N/2-1X2 matrix containing frequency transform of h1 and h2

%% Construct QAM symbol stream.
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');
qamStream = qam_mod(bitStream, M);

% Construct train block.
train_bits = randi([0 1], log2(M)*(N/2-1), 1);
train_frame = qam_mod(train_bits, M);

%% OFDM modulation.
[Tx, a, b, nbPackets] = ofdm_mod_stereo(qamStream, N, Lcp, H, Lt, Ld, train_frame, Equalization);

%% Transmit symbol.
Rx = fftfilt(h1, Tx(:,1)) + fftfilt(h2, Tx(:,2));
Rx = awgn(Rx, SNR, "measured");

%% OFDM demodulation.
[rec_qamStream, CHANNELS] = ofdm_demod_stereo(Rx, N, Lcp, train_frame, Lt, Ld, M, nbPackets, Equalization, mu, alpha);

%% QAM demodulation.
rx_bits = qam_demod(rec_qamStream, M, length(bitStream));

%% Calculate BER
BER = ber(bitStream, rx_bits)