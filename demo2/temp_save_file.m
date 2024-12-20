% DMT-OFDM transmission scheme
%% Cleanup
clear; close all; clc;

%% Parameters
SNRs = 15:5:35;
M = 16; % QAM constellation size
L = 300; % user defined channel order
H = 5; % number of impulse responses
channels = rand(H*(L+1), 1); % Impulse response of channel
N = 2048; % Total number of symbols in a single OFDM frame, i.e., the DFT size
Lcp = max(10, L+1); % Cyclic prefix length, chosen to be longer than the channel impulse response
ON_OFF_mask = ones(N/2-1, 1);

%% Channel effect experiment
% Convert BMP image to bitstream
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');

% Calculate BER for each constellation-SNR combination.
berTransmissions = zeros(H, length(SNRs)); % Placeholder for BERs
for h=1:H % Loop across different transfer functions
    % QAM modulation
    qamStream = qam_mod(bitStream, M);
    
    % QAM constellation visualization
    %scatterplot(qamStream);
    
    % OFDM modulation
    ofdmStream = ofdm_mod(qamStream, N, Lcp, ON_OFF_mask);

    % Channel
    channel = channels((h-1)*L+1:h*L);
    for SNRidx = 1:length(SNRs) % Loop across SNR values
        % Receiving signals
        rxOfdmStream = fftfilt(channel, ofdmStream);
        SNR = SNRs(SNRidx);
        rxOfdmStream = awgn(rxOfdmStream, SNR, "measured");
        
        % OFDM demodulation
        rxQamStream = ofdm_demod(rxOfdmStream, N, Lcp, length(qamStream), channel, ON_OFF_mask, 1);
        
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

%%
N = 2048;
load("channel_session4.mat"); %'acoustic channel impulse response here';
channel = h;
CHANNEL = fft(channel, N);
CHANNEL = CHANNEL(2:N/2);
BWusage = 0.4:0.1:1;
sorted = sort(abs(CHANNEL), 'descend');
Lcp = length(channel);
% Calculate BER for each BWusage.
berTransmissions = zeros(1, length(BWusage)); % Placeholder for BERs
for i=1:length(BWusage)
    % QAM modulation
    qamStream = qam_mod(bitStream, M);

    % ON_OFF_mask
    idx = max(floor(BWusage(i)*length(CHANNEL)), 1);
    threshold = sorted(idx);
    ON_OFF_mask = abs(CHANNEL) >= threshold;

    % OFDM modulation
    ofdmStream = ofdm_mod(qamStream, N, Lcp, ON_OFF_mask);

    % Receiving signals
    rxOfdmStream = fftfilt(channel, ofdmStream);
    SNR = 15;
    rxOfdmStream = awgn(rxOfdmStream, SNR, "measured");

    % OFDM demodulation
    rxQamStream = ofdm_demod(rxOfdmStream, N, Lcp, length(qamStream), channel, ON_OFF_mask, 1);

    % QAM demodulation
    rxBitStream = qam_demod(rxQamStream, M, length(bitStream));

    % Compute BER
    berTransmissions(1,i) = ber(bitStream,rxBitStream);
end