function [OFDM_seq, nbOFDMsymb] = ofdm_mod_pilots(QAM_seq, N, Lcp, trainblock)
% OFDM modulation baed on pilots.
%
% INPUT:
% QAM_seq       T1X1        Sequence containing QAM symbols
% N             1X1         Total number of symbols in a single OFDM frame.
% Lcp           1X1         Cyclic prefix length [samples] (you can ignore this until exercise 3.2.4)
% trainblock    T2X1        Training block of T2 QAM symbols 
%
% OUTPUT:
% data_seq      T3X1        QAM sequence of T3 symbols.
% nbOFDMsymb    1X1         Number of OFDM frames.

%% Extract data
% Padding to make the QAM_seq a multiple of used bins
bins = N/4-1;
padlength = bins - mod(length(QAM_seq), bins);
if padlength ~= bins
    QAM_seq = [QAM_seq; zeros(padlength, 1)];
end

QAM_matrix = zeros(N/2-1, length(QAM_seq)/(bins));
idx = [repmat([0;1], bins, 1); 0];
QAM_matrix(logical(idx)) = reshape(QAM_seq, bins, []); % even bins for data    

% Construct the OFDM frames
fOFDM_frame = [zeros(1, size(QAM_matrix, 2)); QAM_matrix; zeros(1, size(QAM_matrix, 2)); conj(flip(QAM_matrix))];

% Apply the inverse Fourier transform (IFFT)
OFDM_frame = ifft(fOFDM_frame);

% Add in the cyclic prefix
OFDM_frame = [OFDM_frame(end-Lcp+1:end, :); OFDM_frame];
nbOFDMsymb = size(OFDM_frame, 2);

% Serialize the set of OFDM frames
OFDM_seq = reshape(OFDM_frame, [], 1);
end

