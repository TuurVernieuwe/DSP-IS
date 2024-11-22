% Transmission of a picture across the channel
%% Cleanup
clear; close all; clc;

%% Parameters.
N = ; % Total number of symbols in a single OFDM frame, i.e., the DFT size
Lcp = ; % Cyclic prefix length [samples].
M = ; % QAM constellation size.
SNR = ; % SNR of transmission [dB].
Lt = ; % Number of training frames.
Ld = ; % Number of data frames.
fs = ; % Sampling frequency [Hz].
channel = "simulation"; % simulation or acoustic
Equalization = "packet"; % Equalization mode (see ofdm_demod_stereo.m)
Nswitch = (Lt+Ld)*(N+Lcp); % The simulated channel changes every Nswitch number of samples.
smoothing_factor = .99; % Smoothing factor for simulated channel (see simulate_channel.m)

%% Initial channel estimation.
% Construct train block.
train_bits = randi();
train_frame = qam_mod();

% Initial channel estimation.
ofdm_train_seq = ofdm_mod();

% Construct training sequence for transmitter side channel estimation.
ofdm_train_seq_stereo = ;

%% Transmit OFDM sequence.
if channel == "simulation"
    [] = initparams_stereo();
    Rx = simulate_channel_stereo(simin, Nswitch,'channel_stereo_session7.mat',smoothing_factor);
    Rx = awgn(Rx,SNR,'measured');
    aligned_Rx = alignIO();
elseif channel == "acoustic"
    [] = initparams_stereo();
    sim('recplay_stereo');
    Rx = simout.signals.values(:,1);
    aligned_Rx = ;
end

%% Ofdm channel estimation.
CHANNELS = ; % FRequency domain channels
channels = ; % Time domain channels

% Plot channel estimations.
figure(1);
subplot(2,1,1);
plot();
xlabel('Samples')
ylabel('Magnitude [arb.]')   
title('Impulse response.');
legend('left','right');
subplot(2,1,2);
plot();
title('Frequency response.');
xlabel('Samples')
ylabel('Magnitude [dB]')   
legend('left','right');

%% Transmission, including beamforming.
% Construct QAM symbol stream.
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');
qamStream = qam_mod( );

% OFDM modulation.
[Tx, a, b, nbPackets] = ofdm_mod_stereo();

% Transmit data sequence.
if channel == "simulation"
    [] = initparams_stereo();
    Rx = simulate_channel_stereo(simin, Nswitch,'channel_stereo_session7.mat',smoothing_factor);
    Rx = awgn(Rx,SNR,'measured');
    aligned_Rx = alignIO();
    aligned_Rx = aligned_Rx(1:length(Tx)).';
end

%% Post transmission processing
% OFDM demodulation.
[] = ofdm_demod_stereo();
CHANNEL_combo = ; % Extract first estimated channel

% QAM demodulation.
rx_bits = qam_demod();

% Display results Construct image from bitstream
figure;
imageRx = bitstreamtoimage(rx_bits, imageSize, bitsPerPixel);
% Plot images
subplot(1,2,1); colormap(colorMap); image(imageData); axis image; title('Original image'); drawnow;
subplot(1,2,2); colormap(colorMap); image(imageRx); axis image; title(['Received image']); drawnow;
% Calculate BER
BER = ber()