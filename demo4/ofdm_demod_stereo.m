function [ data_seq, CHANNELS] = ofdm_demod_stereo( OFDM_seq, N, Lcp, trainblock, Lt, Ld, M, nbPackets, Equalization, mu, alpha )
% Performs stereo OFDM demodulation
% 
% INPUT:
% OFDM_seq      T1X1            Time domain OFDM sequence of length T samples.
% N             1X1             Total number of symbols in a single OFDM frame.
% Lcp           1X1             Cyclic prefix length [samples].
% trainblock    T2X1            Training block of T2 QAM symbols 
% Lt            1X1             Number of training frames 
% Ld            1X1             Number of data frames 
% M             1X1             QAM-ary constellation size
% Lt            1X1             Number of training frames
% Ld            1X1             Number of data frames
% nbPackets     1X1             Number of packets, where one packet consist of
%                               a training and data frame.
% Equalization  String          Equalization mode: If "fixed" (for transmit_pic_stereo_a), training frames are
%                               transmitted first to estimate the channel and
%                               no additional updating of the channel is
%                               performed thereafter as only data frames
%                               follow. If "packet" packet-based equalization should be
%                               used. (for transmit_pic_stereo_b)
% mu            1X1             NLMS stepsize
% alpha         1X1             NLMS regularization factor
% 
% OUTPUT:
% data_seq      T3X1            QAM sequence of T3 symbols.
% CHANNELS      N/2-1XP         Frequency domain estimated combined channel for each 1..P. 

end

