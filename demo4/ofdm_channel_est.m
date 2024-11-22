function [ CHANNELS ] = ofdm_channel_est( OFDM_seq, N, Lcp, trainblock, Lt )
% Stereo OFDM channel estimation
% 
% INPUT:
% QAM_seq       T1X1        Sequence containing QAM symbols
% N             1X1         Total number of symbols in a single OFDM frame.
% Lcp           1X1         Cyclic prefix length [samples]
% trainblock    T2X1        Training block of T2 QAM symbols
% Lt            1X1         Number of training frames
% 
% OUTPUT:
% CHANNELS      N/2-1X2     Two channel frequency transformed estimated channels

% Reshape the received OFDM sequence (serial to parallel conversion)
OFDM_matrix = ;

% Remove the cyclic prefix
OFDM_matrix = ;

% Apply fft operation
QAM_matrix = ;

% Remove the redundant parts of QAM_matrix
QAM_matrix1 = ;
QAM_matrix2 = ;

% Channel estimation
CHANNELS1 = ;
CHANNELS2 = ;
for n = 1:N/2-1
    CHANNELS1(n) = ;
    CHANNELS2(n) = ;
end

% Concatenate CHANNELS1 and CHANNELS2
CHANNELS = ;
end

