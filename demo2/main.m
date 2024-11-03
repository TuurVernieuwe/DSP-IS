% DMT-OFDM transmission scheme
%% Cleanup
clear; close all; clc;

%% Parameters
M = 16; % QAM constellation size
L = 2; % user defined channel order
channel = rand(1, L); % Impulse response of channel
N = 50; % Total number of symbols in a single OFDM frame, i.e., the DFT size
Lcp = 10; % Cyclic prefix length, chosen to be longer than the channel impulse response
SNR = 30;

%% Channel effect experiment
% Convert BMP image to bitstream
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');

% QAM modulation
qamStream = qam_mod(bitStream, M);

% OFDM modulation
ofdmStream = ofdm_mod(qamStream, N, Lcp);

% Channel
rxOfdmStream = fftfilt(channel, ofdmStream);

% OFDM demodulation
rxQamStream = ofdm_demod(rxOfdmStream, N, Lcp, length(qamStream), channel.');

% QAM demodulation
rxBitStream = qam_demod(rxQamStream, M, length(bitStream));

% Compute BER
% berTransmission = ber(bitStream,rxBitStream);

% Construct image from bitstream
imageRx = bitstreamtoimage(rxBitStream, imageSize, bitsPerPixel);

% Plot images
subplot(2,1,1); colormap(colorMap); image(imageData); axis image; title('Original image'); drawnow;
subplot(2,1,2); colormap(colorMap); image(imageRx); axis image; title(['Received image']); drawnow;
