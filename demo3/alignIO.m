function [ out_aligned ] = alignIO( out, fs )
% Aligns output based on a syncronization pulse.
%
% INPUT:
% out           T1x1    Recorded stream of T1 samples.
% fs            1X1     Sampling frequency [Hz].
% 
% OUTPUT:
% out_aligned   T2X1    Aligned output of T2 samples.

%% Define syncronization pulse
f_sync = 5000;
t = 0:1/fs:1;
sync_pulse = sin(2*pi*f_sync*t)';

%% Align I/O
safety_margin = 50; % Safety margin [samples]

[r, lags] = xcorr(out, sync_pulse); % Find cross-correlation between output and synchronisation pulse
[~,maxIdx] = max(r); % Find index of maximum cross-correlation
delay = lags(maxIdx); % Find delay corresponding to that lag of maximum cross-correlation
startIdx = delay + length(sync_pulse) + fs/10 - safety_margin; % Find start index
out_aligned = out(startIdx:end); % Align output
end

