% Training frames for channel estimation and equalization 
%% Cleanup
clear; clc; close all;

%% Parameters
fs = 16000; % Sampling frequency [Hz]
N = 1024; % Total number of symbols in a single OFDM frame, i.e., the DFT size
Lcp = 600; % Cyclic prefix length [samples]
M = 16; % QAM constellation size
SNR = 35; % SNR of transmission [dB]
accoustic_transmission = 0; % If 1 acoustic transmission occurs, if 0 a simulated transmission.
t = 0:1/fs:1;
f_sync = 5000;
ON_OFF_mask = [ones(N/2-1, 1)];

%% Construct train block.
train_bits = randi([0 1], log2(M)*sum(ON_OFF_mask), 1); % Generate a random vector of bits corresponding to a single OFDM frame
train_block = qam_mod(train_bits, M); % QAM modulate that frame
train_stream = repmat(train_block, 100, 1); % Repeat the training QAM frame
Tx = ofdm_mod(train_stream, N, Lcp, ON_OFF_mask); % OFDM modulate the QAM stream

%% Transmit train block.
if ~accoustic_transmission % Exercise 5.1
    %load('channel_session5.mat');
    load('channel_session4.mat')
    aligned_Rx = fftfilt(h, Tx);
    aligned_Rx = awgn(aligned_Rx, SNR, "measured");
else % Exercise 5.2
    toplay = Tx;
    [simin,nbsecs,fs] = initparams(toplay, fs, sin(2*pi*t*f_sync)');
    sim('recplay');
    Rx = simout.signals.values(:,1);
    aligned_Rx = alignIO(Rx, fs); % Align input and output
    aligned_Rx = aligned_Rx(1:length(Tx));
end

%% OFDM Demodulate
[qam_stream, CHANNELS] = ofdm_demod(aligned_Rx, N, Lcp, length(train_stream), ON_OFF_mask, train_block);

%% QAM Demodulate
rx_bits = qam_demod(qam_stream, M, 100*length(train_bits));

%% BER
BER = ber(repmat(train_bits, 100, 1), rx_bits)

%% Plot (real and) estimated channel.
if ~accoustic_transmission
    figure;
    subplot(2,1,1)
    plot(1:length(h), h);
    title('Real impulse response')
    xlabel('k')
    ylabel('h[k]')
    subplot(2,1,2)
    H = fft(h, N);
    plot(1:N/2-1, pow2db(abs(H(2:N/2)).^2));
    title('Real frequency response')
    xlabel('carrier n')
    ylabel('H(n)')
    xlim([1 N/2-1])

    figure;
    subplot(2,1,1)
    est_h = ifft([0; CHANNELS; 0; flip(conj(CHANNELS))]);
    plot(real(est_h(1:length(h))));
    title('Estimated impulse response')
    xlabel('k')
    ylabel('h[k]')
    subplot(2,1,2)
    plot(1:N/2-1, pow2db(abs(CHANNELS).^2));
    title('Estimated frequency response')
    xlabel('carrier n')
    ylabel('H(n)')
    xlim([1 N/2-1])

    figure;
    subplot(2, 1, 1)
    plot(0:length(h)-1, h-est_h(1:length(h)));
    title('Difference in impulse response')
    xlabel('k')
    ylabel('\Deltah[k]')
    subplot(2, 1, 2)
    plot(1:N/2-1, real(H(2:N/2)) - real(CHANNELS))
    title('Difference in frequency response')
    xlabel('carrier n')
    ylabel('\DeltaH(n)')
    xlim([1 N/2-1])
else
    figure;
    subplot(2,1,1)
    est_h = ifft([0; CHANNELS; 0; flip(conj(CHANNELS))]);
    plot(0:length(est_h)-1, est_h);
    title('Estimated impulse response')
    xlabel('k')
    ylabel('h[k]')
    subplot(2,1,2)
    plot(1:N/2-1, pow2db(abs(CHANNELS).^2));
    title('Estimated frequency response')
    xlabel('carrier n')
    ylabel('H(n)')
end
