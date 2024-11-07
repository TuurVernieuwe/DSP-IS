%% Cleanup
clear; close all; clc;

%% Parameters
M = 16; % QAM constellation size
Nq = log2(M);
h0 = 1; h1 = 0; h2 = 0; % Channel impulse response coefficients
channel = [h0, h1, h2]; % Impulse response of channel
N = 50; % Total number of symbols in a single OFDM frame, i.e., the DFT size
Lcp = 10; % Cyclic prefix length, chosen to be longer than the channel impulse response
SNR = 30;

%% Channel effect experiment
% bit stream 
bitSeq = randi([0, 1], Nq*((N/2)-1), 1); % Generate the bit sequence corresponding to 1 QAM symbol
% bitStream = repmat(bitSeq, N/2-1, 1); % Repeat this bit sequence to fill 1 OFDM frame
bitStream = bitSeq;

% QAM modulation
qamStream = qam_mod(bitStream, M); % Should be of length N/2-1X1. Each symbol should be the same

% QAM constellation visualization
scatterplot(qamStream); % 

% OFDM modulation
ofdmStream = ofdm_mod(qamStream, N, Lcp); % Should be of length N+LcpX1

% Channel
rxOfdmStream = fftfilt(channel, ofdmStream);
rxOfdmStream = awgn(rxOfdmStream, SNR, "measured");

% OFDM demodulation
rxQamStream = ofdm_demod(rxOfdmStream, N, Lcp, length(qamStream));

% QAM constellation visualization
scatterplot(rxQamStream);

% QAM demodulation
rxBitStream = qam_demod(rxQamStream, M, length(bitStream));

% Compute BER
berTransmission = ber(bitStream, rxBitStream);
