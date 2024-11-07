function [ OFDM_seq, nbPackets ] = ofdm_mod( QAM_seq, N, Lcp, varargin)
% OFDM modulation
%
% INPUT:
% QAM_seq       T1X1        Sequence containing QAM symbols
% N             1X1         Total number of symbols in a single OFDM frame.
% Lcp           1X1         Cyclic prefix length [samples] (you can ignore this until exercise 3.2.4)
% varargin
% 
% Session 4 and 5
% ON_OFF_mask   (N/2-1)X1   Mask denoting the bins to use in range (X_1...X_(N/2-1)) with a 1 denothing
%                           that a bin should be used and a 0 denoting that a bin should not be used.
%                           (you can ignore this until exercise 4.3)
%
% Session 6
% ON_OFF_mask   (N/2-1)X1   Mask denoting the bins to use in range (X_1...X_(N/2-1)) with a 1 denothing
%                           that a bin should be used and a 0 denoting that a bin should not be used.
%                           (you can ignore this until exercise 4.3)
% Lt            1X1         Number of training frames (you can ignore this
%                           until exercise 6.1.1)
% Ld            1X1         Number of data frames (you can ignore this
%                           until exercise 6.1.1)
% trainblock    T2X1        Training block of T2 QAM symbols (you can
%                           ignore this until exercise 6.1.1)
% Session 7
% Lt            1X1         Number of training frames 
% trainblock    T2X1        Training block of T2 QAM symbols
% ON_OFF_mask   (N/2-1)X1   Mask denoting the bins to use in range (X_1...X_(N/2-1)) with a 1 denothing
%                           that a bin should be used and a 0 denoting that a bin should not be used.
%
% OUTPUT:
% OFDM_seq      T3X1        Time domain OFDM sequence of length T3 samples.
% nbPackets     1X1         Number of packets, where one packet consist of
%                           a training and data frame. (you can ignore 
%                           this until exercise 6.1.1)

%% Extract data
if nargin == 4
    ON_OFF_mask = varargin{1};
elseif nargin == 7
    ON_OFF_mask = varargin{1};
    Lt = varargin{2};
    Ld = varargin{3};
    trainblock = varargin{4};
elseif nargin == 6
    Lt = varargin{1};
    trainblock = varargin{2};
    ON_OFF_mask = varargin{3};
end

%% Construct the OFDM sequence
% Put the QAM symbols into matrix of N/2-1 rows
bins = sum(ON_OFF_mask);
padlength = bins - mod(length(QAM_seq), bins);
if padlength ~= bins
    QAM_seq = [QAM_seq; zeros(padlength, 1)];
end

QAM_matrix_small = reshape(QAM_seq, bins, []); 
QAM_matrix = zeros(N/2-1, size(QAM_matrix_small, 2));

% Apply on-off mask (you can ignore this until exercise 4.3)
row_small = 1;
for i = 1:length(ON_OFF_mask)
    if ON_OFF_mask(i)
        QAM_matrix(i, :) = QAM_matrix_small(row_small, :);
        row_small = row_small + 1;
    end
end


% Construct the OFDM frames according to Figure 2 in session 3
fOFDM_frame = [zeros(1, size(QAM_matrix, 2)); QAM_matrix; zeros(1, size(QAM_matrix, 2)); conj(flip(QAM_matrix))];

% Apply the inverse Fourier transform (IFFT)
OFDM_frame = ifft(fOFDM_frame);

% Add in the cyclic prefix
OFDM_frame = [OFDM_frame(end-Lcp+1:end, :); OFDM_frame];

% Serialize the set of OFDM frames
OFDM_seq = reshape(OFDM_frame, [], 1);

end
