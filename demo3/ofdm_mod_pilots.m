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
QAM_matrix = ; % even bins for data    

% Construct the OFDM frames
fOFDM_frame = ;

% Apply the inverse Fourier transform (IFFT)
OFDM_frame = ;

% Add in the cyclic prefix
OFDM_frame = ;

% Serialize the set of OFDM frames
OFDM_seq = ;
end

