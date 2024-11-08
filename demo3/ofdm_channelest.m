% Training frames for channel estimation and equalization 
%% Cleanup
clear; clc; close all;

%% Parameters
fs = ; % Sampling frequency [Hz]
N = ; % Total number of symbols in a single OFDM frame, i.e., the DFT size
Lcp = ; % Cyclic prefix length [samples]
M = ; % QAM constellation size
SNR = ; % SNR of transmission [dB]
accoustic_transmission = ; % If 1 acoustic transmission occurs, if 0 a simulated transmission.

%% Construct train block.
train_bits = ; % Generate a random vector of bits corresponding to a single OFDM frame
train_block = qam_mod(); % QAM modulate that frame
train_stream = repmat(); % Repeat the training QAM frame
Tx = ofdm_mod(); % OFDM modulate the QAM stream

%% Transmit train block.
if ~accoustic_transmission % Exercise 5.1
    load('channel_session5.mat');
    aligned_Rx = fftfilt();
    aligned_Rx = awgn();
else % Exercise 5.2
    [] = initparams();
    sim('recplay');
    Rx = simout.signals.values(:,1);
    aligned_Rx = alignIO(); % Align input and output
end

%% OFDM Demodulate
[] = ofdm_demod();

%% QAM Demodulate
rx_bits = qam_demod();

%% BER
BER = ber()

%% Plot (real and) estimated channel.
if ~accoustic_transmission
    figure;
    subplot(2,1,1)
    plot();
    title('Real impulse response')
    xlabel('')
    ylabel('')
    subplot(2,1,2)
    plot();
    title('Real frequency response')
    xlabel('')
    ylabel('')    
end

figure;
subplot(2,1,1)
plot();
title('Estimated impulse response')
xlabel('')
ylabel('')
subplot(2,1,2)
plot();
title('Estimated frequency response')
xlabel('')
ylabel('')