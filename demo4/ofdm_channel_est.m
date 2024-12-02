function [ CHANNELS ] = ofdm_channel_est( OFDM_seq, N, Lcp, trainblock, Lt )
% Stereo OFDM channel estimation
% 
% INPUT:
% OFDM_seq      T1X1        Sequence containing OFDM symbols
% N             1X1         Total number of symbols in a single OFDM frame.
% Lcp           1X1         Cyclic prefix length [samples]
% trainblock    T2X1        Training block of T2 QAM symbols
% Lt            1X1         Number of training frames
% 
% OUTPUT:
% CHANNELS      N/2-1X2     Two channel frequency transformed estimated channels

% Reshape the received OFDM sequence (serial to parallel conversion)
OFDM_matrix = reshape(OFDM_seq, N+Lcp, []);

% Remove the cyclic prefix
OFDM_matrix = OFDM_matrix(Lcp+1:end, :);

% Apply fft operation
QAM_matrix = fft(OFDM_matrix);

% Remove the redundant parts of QAM_matrix
QAM_matrix1 = QAM_matrix(2:size(QAM_matrix, 1)/2, 1:Lt);
QAM_matrix2 = QAM_matrix(2:size(QAM_matrix, 1)/2, Lt+1:end);

% Channel estimation
CHANNELS1 = zeros(N/2-1, 1);
CHANNELS2 = zeros(N/2-1, 1);
for n = 1:N/2-1
    CHANNELS1(n) = (trainblock(n)*ones(Lt, 1)) \ QAM_matrix1(n,:).';
    CHANNELS2(n) = (trainblock(n)*ones(Lt, 1)) \ QAM_matrix2(n,:).';
end

% Concatenate CHANNELS1 and CHANNELS2
CHANNELS = [CHANNELS1, CHANNELS2];
end

