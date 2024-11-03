%% Cleanup
clear; clc; close all;

%% Parameters
SNRs = 1:10; % List of SNR values to consider [dB]
M = 16; % QAM constellation size
N_QAM = 1:6; % Number of symbols per OFDM frame, i.e., the DFT size
N_OFDM = 8;
L = 100000; % Binary sequence length [samples]
Lcp = 2; % Cyclic prefix length [Samples] (you can ignore this until exercise 3.2.4)

%% OFDM experiment
% Calculate BER for each constellation-SNR combination.
BERs = zeros(length(N_QAM), length(SNRs)); % Placeholder for BERs
for n=N_QAM % Loop across values of N
    for SNRidx = 1:length(SNRs)% Loop across SNR values
        % Initialize simulation parameters.
        M = 2^n; % QAM constellation size
        SNR = SNRs(SNRidx); % SNR
        % Generate a pseudo random binary sequence of a user defined length.
        bit_seq = randi([0 1],L,1);
        
        % Modulate bit sequence.
        QAM_seq = qam_mod(bit_seq,M);
        
        % Modulate QAM sequence using OFDM.
        OFDM_seq = ofdm_mod(QAM_seq, N_OFDM, Lcp);
        
        % Add white Gaussian noise.
        rec_OFDM_seq  = awgn(OFDM_seq, SNR);
        
        % Demodulate OFDM sequence.
        [rec_QAM_seq, ~] = ofdm_demod(rec_OFDM_seq, N_OFDM, Lcp, length(QAM_seq));
        
        % Demodulate QAM sequence.
        rec_bit_seq = qam_demod(rec_QAM_seq, M, L);
        
        % Calculate BER.
        BERs(n, SNRidx) = ber(bit_seq, rec_bit_seq);
    end
end

% Plot results.
semilogy(BERs);
ylabel("BER");
xlabel("Bits per QAM symbol");
grid on

% 2.5: formula
% a) log2(M)
% b) N/2-1
% c) fs/(N+L)