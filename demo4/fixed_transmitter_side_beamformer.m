function [a, b] = fixed_transmitter_side_beamformer(H1, H2)
% Takes the two impulse responses as inputs and outputs the optimal
% per-frequency bin scalars a and b
% INPUT:
% H1    N/2-1X1     Frequency response from loudspeaker 1
% H2    N/2-1X1     Frequency response from loudspeaker 2
%
% OUPUT:
% a     N/2-1X1     Optimal per-frequency bin scalars for loudspeaker 1
% b     N/2-1X1     Optimal per-frequency bin scalars for loudspeaker 2

denom = sqrt(H1.*conj(H1) + H2.*conj(H2));
a = conj(H1)./denom;
b = conj(H2)./denom;

% Plot
figure;
plot([pow2db(abs(H1).^2), pow2db(abs(H2).^2), pow2db(abs(denom).^2)])
axis tight
xlabel('k')
ylabel('Transfer function H(k)')
legend('H^1', 'H^2', 'H^{1+2}')
end