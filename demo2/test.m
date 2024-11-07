% DMT-OFDM transmission scheme
%% Cleanup
clear; close all; clc;

%% Parameters
SNR= 15;
M = 16; % QAM constellation size
load("channel.mat");
L = 300; % user defined channel order
channel = h; % Impulse response of channel
N = 1024; % Total number of symbols in a single OFDM frame, i.e., the DFT size
Lcp = max(10, length(channel)); % Cyclic prefix length, chosen to be longer than the channel impulse response
BWusage = 0.6;
Gamma = 10;

%% Channel effect experiment
% Convert BMP image to bitstream
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');

% Calculate the channel frequency response
H = fft(channel, N);
CHANNEL = H(2:N/2);
sorted = sort(abs(CHANNEL), 'descend');

% ON_OFF_mask
idx = max(floor(BWusage*length(CHANNEL)), 1);
threshold = sorted(idx);
ON_OFF_mask = abs(CHANNEL) >= threshold;

noisePower = 0;
for i = 1:N/2-1
    noisePower = noisePower + abs(CHANNEL(i)).^2;
end

% Calculate bit allocation (adaptive bitloading)
M_k = zeros(length(CHANNEL), 1);
for k = 1:length(CHANNEL)
    b_k = floor(log2(1 + ((N/2-1) * db2pow(SNR) * abs(CHANNEL(k)).^2) / (Gamma * noisePower)));
    if b_k<1
        b_k = 1;
    end
    M_k(k) = 2^b_k;
end

sm = 0;
for k = 1:length(M_k)
    sm = sm + log2(M_k(k));
end

padLength = sm - mod(length(bitStream),sm);
if padLength ~= sm
    bitStream_padded = [bitStream; zeros(padLength, 1)];
end
l = length(bitStream_padded)/sm;

% QAM modulation
qamStream = adaptive_qam_mod(bitStream_padded, M_k, l);

% QAM constellation visualization
scatterplot(qamStream);

% OFDM modulation
ofdmStream = ofdm_mod(qamStream, N, Lcp, ON_OFF_mask);

% Receiving signals
rxOfdmStream = fftfilt(channel, ofdmStream);
rxOfdmStream = awgn(rxOfdmStream, SNR, "measured");

% OFDM demodulation
rxQamStream = ofdm_demod(rxOfdmStream, N, Lcp, length(qamStream), channel, ON_OFF_mask, 1);

% QAM constellation visualization
scatterplot(rxQamStream);

% QAM demodulation
rxBitStream = adaptive_qam_demod(rxQamStream, M_k, length(bitStream), l);

% Compute BER
berTransmission = ber(bitStream,rxBitStream);

% Construct image from bitstream
imageRx = bitstreamtoimage(rxBitStream, imageSize, bitsPerPixel);

% Plot images
figure
subplot(2,1,1); colormap(colorMap); image(imageData); axis image; title('Original image'); drawnow;
subplot(2,1,2); colormap(colorMap); image(imageRx); axis image; title('Received image'); drawnow;