function [OFDM_seq, a, b, nbPackets] = ofdm_mod_stereo( QAM_seq, N, Lcp, H, Lt, Ld, trainblock, Equalization)
% Stereo OFDM modulation
% 
% INPUT:
% QAM_seq       T1X1        Sequence containing QAM symbols
% N             1X1         Total number of symbols in a single OFDM frame.
% Lcp           1X1         Cyclic prefix length [samples].
% H             N/2-1X2     Two frequency domain channels for a DFT size of N.
% Lt            1X1         Number of training frames.
% Ld            1X1         Number of data frames.
% trainblock    T2X1        Training block of T2 QAM symbols
% Equalization  String      Equalization mode: If "fixed" (for transmit_pic_stereo_a), training frames are
%                           transmitted first to estimate the channel and
%                           no additional updating of the channel is
%                           performed thereafter as only data frames
%                           follow. If "packet" packet-based equalization should be
%                           used. (for transmit_pic_stereo_b)
% 
% OUTPUT:
% OFDM_seq      T3X2        Two channel time domain OFDM sequence of length T3 samples.
% a             N/2-1X1     Per-bin factor associated with the first stream.
% b             N/2-1X1     Per-bin factor associated with the second stream.
% nbPackets     1X1         Number of packets, where one packet consist of
%                           a training and data frame.
end