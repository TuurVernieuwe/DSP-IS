function [ bit_seq ] = adaptive_qam_demod( QAM_seq, M_values, streamLength, frame_size)
% Demodulates M-ary QAM symbols to bits. 
%
% INPUT:
% QAM_seq       T2X1    Modulated bit sequence into M-aray QAM format of length
%                       T2 symbols.
% M             1X1     M-aray QAM format (corresponding to an integer power of 2)
% varargin:
% streamLength  1X1     Length of returned bit sequence [samples].
%
% OUTPUT:
% bit_seq       T1x1    Bit sequence of T1 bits 

bit_seq = [];
index = 1;

for k = 1:length(M_values)
    M = M_values(k);
    N = log2(M);
    QAMs_k = QAM_seq(index:index+frame_size-1);
    bits = qamdemod(QAMs_k, M, 'OutputType', 'bit', UnitAveragePower=true);
    bit_seq = [bit_seq; bits];
    index = index + frame_size;
end

bit_seq = bit_seq(1:streamLength);

end

