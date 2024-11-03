%% Cleanup
clear; clc; close all;

%% Parameters
L = 10507; % Scalar length of binary sequence [samples]
SNRs = 1:10; % List of SNR values to consider [dB]
N = 1:6; % List of number of bits per QAM symbol to consider

%% Calculate BER for each constellation-SNR combination.
BERs = zeros(length(N), length(SNRs)); % Placeholder for BERs
for n=N % Loop across values of N
    for SNRidx = 1:length(SNRs)% Loop across SNR values
        % Initialize simulation parameters.
        M = 2^n; % QAM constellation size
        SNR = SNRs(SNRidx); % SNR
        
        % Generate a pseudo random binary sequence of a user defined length.
        bit_seq = randi([0 M-1],L,1);
        
        % Modulate bit sequence.
        QAM_seq = qam_mod(bit_seq,M);
        
        % Add white Gaussian noise.
        req_QAM_seq  = awgn(QAM_seq, SNR);        

        % Demodulate QAM sequence.
        rec_bit_seq = qam_demod(req_QAM_seq, M, L);
        
        % Calculate BER.
        BERs(n, SNRidx) = ber(bit_seq, rec_bit_seq);
    end
end

% Plot results.
semilogy(BERs);
ylabel("BER");
xlabel("Bits per QAM symbol");
grid on