function [ BER ] = ber( input1,input2 )
% Computes the bit error rate (BER) between two bit sequences input1 and
% input2.
%
% INPUT:
% input1    TX1     Bit sequence 1 of T bits.
% input2    TX1     Bit sequence 2 of T bits.
% 
% OUTPUT:
% BER       1X1     Bit error rate, in [0,1].

%% Compute the BER
[~, BER] = biterr(input1, input2);

%% Assert that the BER is a scalar value
assert(isscalar(BER),'BER is not a scalar.')

end
