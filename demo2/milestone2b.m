% DMT-OFDM transmission scheme
%% Cleanup
clear; close all; clc;

%% Parameters
SNR= 15;
M = 16; % QAM constellation size
load("channel.mat");
L = 300; % user defined channel order
channel = h; % Impulse response of channel
N = 2048; % Total number of symbols in a single OFDM frame, i.e., the DFT size
Lcp = max(10, length(channel)); % Cyclic prefix length, chosen to be longer than the channel impulse response
BWusage = 0.6;

%% Channel effect experiment
% Convert BMP image to bitstream
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');

CHANNEL = fft(channel, N);
CHANNEL = CHANNEL(2:N/2);
sorted = sort(abs(CHANNEL), 'descend');

% ON_OFF_mask
idx = max(floor(BWusage*length(CHANNEL)), 1);
threshold = sorted(idx);
ON_OFF_mask = abs(CHANNEL) >= threshold;

% QAM modulation
qamStream = qam_mod(bitStream, M);

% QAM constellation visualization
%scatterplot(qamStream);

% OFDM modulation
ofdmStream = ofdm_mod(qamStream, N, Lcp, ON_OFF_mask);

% Receiving signals
rxOfdmStream = fftfilt(channel, ofdmStream);
rxOfdmStream = awgn(rxOfdmStream, SNR, "measured");

% OFDM demodulation
rxQamStream = ofdm_demod(rxOfdmStream, N, Lcp, length(qamStream), channel, ON_OFF_mask, 1);

% QAM constellation visualization
%scatterplot(rxQamStream);

% QAM demodulation
rxBitStream = qam_demod(rxQamStream, M, length(bitStream));

% Compute BER
berTransmission = ber(bitStream,rxBitStream);

% Construct image from bitstream
imageRx = bitstreamtoimage(rxBitStream, imageSize, bitsPerPixel);

% Plot images
figure
subplot(2,1,1); colormap(colorMap); image(imageData); axis image; title('Original image'); drawnow;
subplot(2,1,2); colormap(colorMap); image(imageRx); axis image; title('Received image'); drawnow;