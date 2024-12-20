function [data_seq, CHANNELS] = ofdm_demod_pilots( OFDM_seq, N, Lcp, streamLength, trainblock,nbOFDMsymb)
% OFDM demodulation using pilot tones
% INPUT:
% OFDM_seq      T1X1            Time domain OFDM sequence of length T1 samples.
% N             1X1             Total number of symbols in a single OFDM frame.
% Lcp           1X1             Cyclic prefix length [samples] (you can ignore this until exercise 3.2.4)
% streamLength  1X1             Length of QAM sequence after call to qammod [samples].
% varargin  
% 
% Session 5
% trainblock    T2X1            Training block of T2 QAM symbols (you can
%                               ignore this until exercise 6.1.1)
% nbOFDMsymb    1X1             Number of OFDM frames.
%
% OUTPUT:
% data_seq      T3X1            QAM sequence of T3 symbols.
% CHANNELS      N/2-1XnbOFDMsymbFrequency domain estimated channel for each frame. (
%                               you can ignore this until exercise 5.1.4)   

%% Perform OFDM demodulation
% Reshape the received OFDM sequence (serial to parallel conversion)
OFDM_matrix = reshape(OFDM_seq, N+Lcp, []);

% Remove the cyclic prefix (you can ignore this until exercise 3.2.4)
OFDM_matrix = OFDM_matrix(Lcp+1:end, :);

% Apply fft operation
QAM_matrix = fft(OFDM_matrix);

% Remove the redundant parts of QAM_matrix
QAM_matrix = QAM_matrix(2:size(QAM_matrix, 1)/2, :);

data_matrix = zeros(N/4-1, size(QAM_matrix, 2)); % Placeholder for data
CHANNELS = zeros(N/4, 1); % Placeholder for channels

for pIdx = 1:nbOFDMsymb % Loop across frames
    % Extract packet.
    idx = [repmat([0;1], bins, 1); 0];
    dataFrames(:,pIdx) = QAM_matrix(logical(idx),pIdx);
        
    % Channel estimation
    CHANNELS(:,pIdx) = QAM_matrix(~logical(idx),pIdx); % Save channel
    
    % Equalization of data frames
    data_matrix(:,pIdx) = ;
end

data_seq = ; % Return data sequence





