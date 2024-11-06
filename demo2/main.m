% DMT-OFDM transmission scheme
%% Cleanup
clear; close all; clc;

%% Parameters
SNRs = 10:5:30;
M = 16; % QAM constellation size
L = 2; % user defined channel order
H = 5; % number of impulse responses
channels = rand(1, H*(L+1)); % Impulse response of channel
N = 50; % Total number of symbols in a single OFDM frame, i.e., the DFT size
Lcp = max(10, L+1); % Cyclic prefix length, chosen to be longer than the channel impulse response

%% Channel effect experiment
% Convert BMP image to bitstream
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');

% Calculate BER for each constellation-SNR combination.
berTransmissions = zeros(H, length(SNRs)); % Placeholder for BERs
for h=1:H % Loop across different transfer functions
    for SNRidx = 1:length(SNRs)% Loop across SNR values
        % QAM modulation
        qamStream = qam_mod(bitStream, M);
        
        % QAM constellation visualization
        %scatterplot(qamStream);
        
        % OFDM modulation
        ofdmStream = ofdm_mod(qamStream, N, Lcp);
        
        % Channel
        channel = channels((h-1)*L+1:h*L);
        rxOfdmStream = fftfilt(channel, ofdmStream);
        SNR = SNRs(SNRidx);
        rxOfdmStream = awgn(rxOfdmStream, SNR);
        
        % OFDM demodulation
        rxQamStream = ofdm_demod(rxOfdmStream, N, Lcp, length(qamStream), channel.');
        
        % QAM constellation visualization
        %scatterplot(rxQamStream);
        
        % QAM demodulation
        rxBitStream = qam_demod(rxQamStream, M, length(bitStream));
        
        % Compute BER
        berTransmissions(h,SNRidx) = ber(bitStream,rxBitStream);
    end
end

% Plot results.
plot(1:H, berTransmissions);
ylabel("BER");
xlabel("Transfer function H(z)");
set(gca,'xtick',1:H)
title('SNRs = 10:5:30')
grid on

% Construct image from bitstream
imageRx = bitstreamtoimage(rxBitStream, imageSize, bitsPerPixel);

% Plot images
figure
subplot(2,1,1); colormap(colorMap); image(imageData); axis image; title('Original image'); drawnow;
subplot(2,1,2); colormap(colorMap); image(imageRx); axis image; title('Received image'); drawnow;
