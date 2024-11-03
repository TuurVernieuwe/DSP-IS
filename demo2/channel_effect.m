%% Cleanup
clear; close all; clc;

%% Parameters
M = ; % QAM constellation size
channel = ; % Impulse response of channel
N = ; % Total number of symbols in a single OFDM frame, i.e., the DFT size

%% Channel effect experiment
% bit stream 
bitSeq = ; % Generate the bit sequence corresponding to 1 QAM symbol
bitStream = repmat(); % Repeat this bit sequence to fill 1 OFDM frame

% QAM modulation
qamStream = qam_mod; % Should be of length N/2-1X1. Each symbol should be the same

% QAM constellation visualization
scatterplot(); % 

% OFDM modulation
ofdmStream = ofdm_mod( ); % Should be of length N+LcpX1

% Channel
rxOfdmStream = fftfilt();
rxOfdmStream = awgn();

% OFDM demodulation
rxQamStream = ofdm_demod( );

% QAM constellation visualization
scatterplot();

% QAM demodulation
rxBitStream = qam_demod( );

% Compute BER
berTransmission = ber();
