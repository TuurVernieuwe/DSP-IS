% Decision direction mode of adaptive linear filtering
%% Cleanup
clc; clear; close all;

%% Parameters.
M = 16; % QAM constelllation size
nbQAMsymb = 1000; % Number of QAM symbols 
Hk = randn(1,1,'like',1i); % Channel to consider (Only one frequency bin will be considered here)
% Hk = 0.8 + 0.6j;
alpha = 1; % Regularisation constant
SNR = 30; % Signal-to-noise-ratio [dB]

%% Construct data sequences.
nbBits = log2(M)*nbQAMsymb ; % Number of bits required to generate nbQAMsymb symbols
bitstream = randi([0 1], nbBits, 1); % bitstream of nbBits bits
Xk = qam_mod(bitstream, M); % QAM symbol sequence
Yk = awgn(Hk*Xk, SNR); % Recorded QAM symbol sequence
mu = [0.02, 0.1, 0.2, 0.5, 1, 5, 10]; % List of stepsizes
legendCell = cell(1, length(mu)); % Preallocate cell array

for iN=1:length(mu) % Counter of stepsizes
    % NMLS filter implementation.
    % Initialise filters, reconstructed transmitted signal and error
    w = zeros(nbQAMsymb, 1); 
    rec_Xk = zeros(nbQAMsymb,1); Ek = zeros(nbQAMsymb,1);
    w(1) = (1 + 0.2)/conj(Hk); 
    for n = 1:nbQAMsymb-1
        % Apply filter.
        estXk = conj(w(n))*Yk(n+1);
        % Reconstruct transmitted signal.
        rec_bits = qam_demod(estXk, M, log2(M));
        rec_Xk(n) = qam_mod(rec_bits, M);
        % Calculate error signal.
        Ek(n) = rec_Xk(n) - estXk;
        % Update filter.
        w(n+1) = w(n) + mu(iN)/(alpha + conj(Yk(n+1))*Yk(n+1)) * Yk(n+1) * conj(Ek(n));
    end
    legendCell{iN} = num2str(mu(iN),'mu=% .2f');
    % Plot results.
    figure(1);
    semilogy(abs(1./w'-Hk));
    title('NMLS estimation error.');
    xlabel('Iteration'); ylabel('Error magnitude');
    hold on
    legend(legendCell{1:iN})
end


