% Parameters
Fs = 16000;               % Sampling frequency (Hz)
duration = 2;            % Duration of the signal (seconds)
N = 512;                % Number of points for the FFT
f_stop1 = 700;          % Lower stop band frequency (Hz)
f_stop2 = 3000;         % Upper stop band frequency (Hz)

% Generate white noise
t = 0:1/Fs:duration;    % Time vector
white_noise = wgn(duration*Fs, 1, 1); % Generate white noise

% Design the band-stop filter
Wn = [f_stop1 f_stop2] / (Fs / 2); % Normalize frequency
b = fir1(100, Wn, 'stop'); % Design the filter using fir1
[H, f] = freqz(b, 1, N, Fs); % Frequency response of the filter

% Filter the white noise
filtered_noise = filter(b, 1, white_noise);

% Plotting
figure;

% Plot white noise
subplot(2, 1, 1);
plot(t(1:1000), white_noise(1:1000)); % Plot the first 1000 samples
title('White Noise Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Plot filtered noise
subplot(2, 1, 2);
plot(t(1:1000), filtered_noise(1:1000)); % Plot the first 1000 samples
title('Band-Stop Filtered Noise Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Frequency response plot
figure;
plot(f, 20*log10(abs(H))); % Convert to dB
title('Frequency Response of the Band-Stop Filter');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;
axis([0 Fs/2 -60 5]); % Set axis limits
