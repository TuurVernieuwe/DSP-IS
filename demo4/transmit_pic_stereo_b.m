% Transmission of a picture across the channel
%% Cleanup
clear; close all; clc;

%% Parameters.
N = 1024; % Total number of symbols in a single OFDM frame, i.e., the DFT size
Lcp = 300; % Cyclic prefix length [samples].
M = 8; % QAM constellation size.
SNR = 15; % SNR of transmission [dB].
Lt = 10; % Number of training frames.
Ld = 10; % Number of data frames.
fs = 16000; % Sampling frequency [Hz].
channel = "simulation"; % simulation or acoustic
Equalization = "packet"; % Equalization mode (see ofdm_demod_stereo.m)
Nswitch = (Lt+Ld)*(N+Lcp); % The simulated channel changes every Nswitch number of samples.
smoothing_factor = .99; % Smoothing factor for simulated channel (see simulate_channel.m)

%% Initial channel estimation.
% Construct train block.
train_bits = randi([0 1], log2(M)*(N/2-1), 1);
train_frame = qam_mod(train_bits, M);

% Initial channel estimation.
ofdm_train_seq = ofdm_mod(repmat(train_frame, Lt, 1), N, Lcp, ones(N/2-1,1));

% Construct training sequence for transmitter side channel estimation.
ofdm_train_seq_stereo = [ofdm_train_seq, zeros(length(ofdm_train_seq), 1); zeros(length(ofdm_train_seq), 1), ofdm_train_seq];

%% Transmit OFDM sequence.
if channel == "simulation"
    [simin, nbsecs, fs] = initparams_stereo(ofdm_train_seq_stereo, fs);
    Rx = simulate_channel_stereo(simin, Nswitch,'channel_stereo_session7.mat',smoothing_factor);
    Rx = awgn(Rx,SNR,'measured');
    aligned_Rx = alignIO(Rx, fs);
    aligned_Rx = aligned_Rx(1:length(ofdm_train_seq_stereo));
elseif channel == "acoustic"
    [simin, nbsecs, fs] = initparams_stereo(ofdm_train_seq_stereo, fs);
    sim('recplay_stereo');
    Rx = simout.signals.values(:,1);
    aligned_Rx = alignIO(Rx, fs);
    aligned_Rx = aligned_Rx(1:length(ofdm_train_seq));
end

%% Ofdm channel estimation.
CHANNELS = ofdm_channel_est(aligned_Rx, N, Lcp, train_frame, Lt); % Frequency domain channels
channels = ifft(CHANNELS); % Time domain channels

% Plot channel estimations.
figure(1);
subplot(2,1,1);
plot(abs(channels));
xlabel('Samples')
ylabel('Magnitude [arb.]')   
title('Impulse response.');
legend('left','right');
subplot(2,1,2);
plot(pow2db(abs(CHANNELS).^2));
xlim([1 N/2-1])
ylim([-60 10])
title('Frequency response.');
xlabel('Samples')
ylabel('Magnitude [dB]')   
legend('left','right');

%% Transmission, including beamforming.
% Construct QAM symbol stream.
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');
qamStream = qam_mod(bitStream, M);

% OFDM modulation.
[Tx, a, b, nbPackets] = ofdm_mod_stereo(qamStream, N, Lcp, CHANNELS, Lt, Ld, train_frame, Equalization);

% Transmit data sequence.
if channel == "simulation"
    [simin, nbsecs, fs] = initparams_stereo(Tx, fs);
    Rx = simulate_channel_stereo(simin, Nswitch,'channel_stereo_session7.mat',smoothing_factor);
    Rx = awgn(Rx,SNR,'measured');
    aligned_Rx = alignIO(Rx, fs);
    aligned_Rx = aligned_Rx(1:length(Tx)).';
end

%% Post transmission processing
% OFDM demodulation.
mu = 0.5; % NLMS stepsize
alpha = 1; % NLMS regularization
[rx_qam, CHANNELS] = ofdm_demod_stereo(aligned_Rx, N, Lcp, train_frame, Lt, Ld, M, nbPackets, Equalization, mu, alpha);
CHANNEL_combo = CHANNELS(:,1); % Extract first estimated channel

% QAM demodulation.
rx_bits = qam_demod(rx_qam, M, length(bitStream));

% Display results Construct image from bitstream
figure;
imageRx = bitstreamtoimage(rx_bits, imageSize, bitsPerPixel);
% Plot images
subplot(1,2,1); colormap(colorMap); image(imageData); axis image; title('Original image'); drawnow;
subplot(1,2,2); colormap(colorMap); image(imageRx); axis image; title('Received image'); drawnow;
% Calculate BER
BER = ber(bitStream, rx_bits)