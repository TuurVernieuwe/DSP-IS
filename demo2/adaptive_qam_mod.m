function [ QAM_seq ] = adaptive_qam_mod( bit_seq, M_values, l)
% Modulates a bit sequence into M-ary QAM format. 
%
% INPUT:
% bit_seq   T1x1    Bit sequence of T1 bits 
% M         1X1     M-aray QAM format (corresponding to an integer power of 2)
%
% OUTPUT:
% QAM_seq   T2X1    Modulated bit sequence into M-aray QAM format of length
%                   T2 symbols.
QAM_seq = [];
N_bits = 0;

for k = 1:length(M_values)
    M = M_values(k);
    assert(sum(nextpow2(M)==log2(M))==length(M),'M is not a power of 2.')
    N = log2(M); % Number of bits per QAM symbol
    

    %Get the bits for this QAM symbol
    bits_k = bit_seq(N_bits + 1: N_bits + N*l);
    N_bits = N_bits + N*l;

    % Call to qammod() to obtain the M-ary QAM symbols
    QAM_symbols = qammod(bits_k, M, 'InputType', 'bit', UnitAveragePower=true);
    QAM_seq = [QAM_seq; QAM_symbols];
end

end

